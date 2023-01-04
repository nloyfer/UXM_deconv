
import os
import sys
import time
import pandas as pd
import os.path as op
import argparse
import subprocess
import multiprocessing
from multiprocessing import Pool
import hashlib
import gzip
import tempfile
from pathlib import Path


dpath = str(Path(op.realpath(__file__)).parent.parent)
DEF_TMP_DIR = op.join(dpath, 'tmp_dir')
DEF_ATLAS = op.join(dpath, 'supplemental/Atlas.U25.l4.hg19.tsv')
# DEF_NR_THREADS = multiprocessing.cpu_count()
DEF_NR_THREADS = 4
coord_cols = ['chr', 'start', 'end', 'startCpG', 'endCpG']

##############################
#                            #
#   General methods          #
#                            #
##############################

def eprint(*args,  **kargs):
    print(*args, file=sys.stderr, **kargs)


def validate_file(fpath):
    if not op.isfile(fpath):
        eprint('Invalid file', fpath)
        exit()
    return fpath


def check_executable(cmd):
    for p in os.environ['PATH'].split(":"):
        if os.access(op.join(p, cmd), os.X_OK):
            return
    eprint(f'executable {cmd} not found in PATH')
    exit(1)


def pat2name(pat):
    return op.basename(op.splitext(op.splitext(pat)[0])[0])


def drop_dup_keep_order(lst):
    seen = set()
    return [x for x in lst if not (x in seen or seen.add(x))]

def remove_files(files):
    for f in files:
        if op.isfile(f):
            os.remove(f)

def mkdir_p(dirpath):
    if not op.isdir(dirpath):
        os.mkdir(dirpath)
    return dirpath


##############################
#                            #
#     Homog logic            #
#                            #
##############################


def clear_mem_file(tmp_dir_l):
    try:
        for d in os.listdir(tmp_dir_l):
            if not d.endswith('.bed'):
                continue
            dd = op.join(tmp_dir_l, d)
            if not op.isfile(dd):
                continue
            if (time.time() - op.getmtime(dd)) > 5 * 60 * 60:
                os.remove(dd)

    except Exception as e:
        # no reason to crash over this cleanup process
        return



def load_homog(homog_path):
    validate_file(homog_path)
    names = coord_cols + list('UXM')
    ddtype = {x: int for x in names}
    ddtype['chr'] = str
    try:
        df = pd.read_csv(homog_path, sep='\t', dtype=ddtype, names=names, comment='#')
    except:
        eprint(f'Failed loading memoization file {homog_path} .')
        # os.remove(homog_path)
        return pd.DataFrame()
    return df


def pat2memfile(pat, tmp_dir_l):
    # get memoization path for a given pat path
    dir_hash = hashlib.md5(op.dirname(op.realpath(pat))\
            .encode()).hexdigest()[:6]
    return op.join(tmp_dir_l, f'{pat2name(pat)}.mem.{dir_hash}.homog.gz')


def pat2homog_mp_wrap(markers, pats, tmp_dir_l, rlen, verb, force, nodump, debug, threads):
    # multiprocess wrapper for the pat2homog method
    # return a dict {pat: homog table}
    # homog table is markers x [U, X, M]
    params = [(markers, p, tmp_dir_l, rlen, verb, force, nodump, debug)
               for p in pats]
    p = Pool(threads)
    arr = p.starmap(pat2homog, params)
    p.close()
    p.join()

    res = {k: v for d in arr for k, v in d.items()}
    # validate results:
    for k in res.keys():
        if res[k].empty:
            eprint('Failed in homog_mem: ', k)
            exit(1)
    return res


def gen_homogs(markers, pats, tmp_dir, verb, rlen, force, nodump, debug, threads):

    check_executable('wgbstools')

    # arguments cleanup
    markers = markers[coord_cols]          # keep only coords columns
    # change to abs path (for mem files), validate existance, drop dups
    pats = sorted(set(map(validate_file, map(op.abspath, pats))))
    # create tmp_dir if needed
    tmp_dir_l = mkdir_p(op.join(mkdir_p(tmp_dir), f'l{rlen}'))

    # run compute homog values on missing markers and pats
    uxm_dict = pat2homog_mp_wrap(markers, pats, tmp_dir_l, rlen, verb, force, nodump, debug, threads)

    # clean tmp_dir_l from old temp files:
    clear_mem_file(tmp_dir_l)
    return uxm_dict


def wrap_cpp_tool(pat, markers, tmp_dir_l, rlen, verb, debug):
    # dump extended markers file (for tabix -R)
    validate_file(pat)

    # dump reducerd markers file
    uniq_name = pat2name(pat) + '.' + next(tempfile._get_candidate_names())
    uniq_name = op.join(tmp_dir_l, uniq_name)
    tmp_mpath = uniq_name + '.bed'
    markers[coord_cols].to_csv(tmp_mpath, sep='\t', header=None, index=None)

    # pat to homog.gz:
    cmd = f'wgbstools homog -f --rlen {rlen} -b {tmp_mpath} {pat} --prefix {uniq_name}'
    if debug:
        cmd += ' -v '
        eprint(cmd)
    so = None if verb else subprocess.PIPE
    subprocess.check_call(cmd, shell=True, stderr=so, stdout=so)
    if not debug:
        remove_files([tmp_mpath])
    return uniq_name + '.uxm.bed.gz'


def pat2homog(markers, pat, tmp_dir_l, rlen, verb, force, nodump, debug):

    mempat = pat2memfile(pat, tmp_dir_l)
    name = pat2name(pat)

    # load current markers if exist:

    omark = markers.copy()

    markers = markers.drop_duplicates(subset=coord_cols).reset_index(drop=True)
    # in case this pat file is unseen before, 
    # or it's mempat is older than the pat, or --force is specified:
    ignore_mem = False
    if not op.isfile(mempat):
        msg = f'[ {name} ] no memoization found'
        ignore_mem = True
    elif op.getctime(mempat) < op.getctime(pat):
        msg = f'[ {name} ] memoization is older than pat! deleting it'
        ignore_mem = True
    if force:
        msg = f'[ {name} ] overwriting existing memoization file (--force)'
        ignore_mem = True
    if ignore_mem:
        remove_files([mempat])
        if verb:
            eprint(msg)
        remain_mrk = markers
        homog = pd.DataFrame()

    # else, find markers not already present in mempat
    else:
        homog = load_homog(mempat).reset_index(drop=True)
        if homog.empty:
            if verb:
                eprint(f'loaded empty homog for {pat}: {mempat}')
            os.remove(mempat)
            return {pat: homog}
        remain_mrk = markers.merge(homog, how='left')
        remain_mrk = remain_mrk[remain_mrk['U'].isna()]
        if verb and not remain_mrk.empty:
            eprint(f'[ {name} ] found memoization, missing {remain_mrk.shape[0]} markers' )

    # if all markers are present, return them
    if remain_mrk.empty:
        if verb:
            eprint(f'[ {name} ] all markers found in memory')
        res = omark.merge(homog, how='left')
        # debug \ validation
        if res.shape[0] != omark.shape[0]:
            eprint(f'[ {name} ] Error:', res.shape, omark.shape)
            return {pat: pd.DataFrame()}
        return {pat: res[['U', 'X', 'M']]}

    # otherwise, run compute remaining homog values 
    tmp_homog_path = wrap_cpp_tool(pat, remain_mrk, tmp_dir_l, rlen, verb, debug)

    # homog.gz to numpy array uxm:
    uxm = load_homog(tmp_homog_path)
    # cleanup 
    remove_files([tmp_homog_path])
    nodump = bool(nodump)
    if uxm[['U', 'X', 'M']].values.sum() == 0:
        if verb:
            eprint('\033[91m' + f'WARNING:' + '\033[0m' +
                    f' possibly failed in {pat} - all {uxm.shape[0]} ' \
                    'values are zero. memoization is not updated, to be safe')
        nodump = True

    all_markers = pd.concat([homog, uxm])
    all_markers.sort_values(by=['startCpG', 'endCpG'], inplace=True)
    res = omark.merge(all_markers, how='left')
    if not nodump:
        # dump
        with gzip.open(mempat, 'wt') as mp:
            mp.write(f'# {pat}\n')
        all_markers.to_csv(mempat, sep='\t', index=None, na_rep='NA', header=None, mode='a')
        # chmod to -rw-rw-r-- so all group users can share the same temp_dir
        os.chmod(mempat, 0o664)

    return {pat: res[['U', 'X', 'M']]}


def deduce_l_from_name(args):
    if args.rlen > 0:
        return args.rlen
    name = op.splitext(op.basename(args.atlas))[0]
    for i in range(3, 6):
        if f'l{i}' in name:
            if args.verbose:
                eprint(f'infered from atlas name that rlen={i}')
            return i
    eprint(f'Error: could not deduce -l value from atlas name {name}')
    exit(1)


def add_memoiz_args(parser):
    parser.add_argument('--tmp_dir', '-T', default=DEF_TMP_DIR,
            help=f'path to directory to store temporary memoization ' \
            f'files [{DEF_TMP_DIR}]')
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--nodump', action='store_true',
            help='Do not update memoization files. This flag is important ' \
                 'when running uxm on multiple machines that share the same' \
                 ' temp_dir')
    parser.add_argument('--force', '-f', action='store_true',
            help='Force run homog: recreate the whole memoization table ')
    parser.add_argument('--threads', '-@', type=int, default=DEF_NR_THREADS,
            help=f'Number of threads [DEF_NR_THREADS]')
    parser.add_argument('--rlen', '-l', type=int, default=-1,
            help='minimal CpGs per read required to consider the read.' \
                    ' By default this value is deduced from the atlas name (e.g. atlas.l4.tsv)')
    parser.add_argument('--debug', '-d', action='store_true')

