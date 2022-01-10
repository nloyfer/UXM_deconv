
import os
import sys
import pandas as pd
import numpy as np
import os.path as op
import argparse
import subprocess
from multiprocessing import Pool
import uuid
import hashlib
import gzip


DEF_TMP_DIR = './tmp_dir'
coord_cols = ['chr', 'start', 'end', 'startCpG', 'endCpG']
cc = ['startCpG', 'endCpG']

##############################
#                            #
#   General methods          #
#                            #
##############################

def eprint(*args,  **kargs):
    print(*args, file=sys.stderr, **kargs)

def validate_file(fpath):
    # if 'Megakar' in fpath: # todo: MK hack
        # return
    if not op.isfile(fpath):
        eprint('Invalid file', fpath)
        exit()

def pat2name(pat):
    return op.basename(op.splitext(op.splitext(pat)[0])[0])

def remove_files(files):
    for f in files:
        if op.isfile(f):
            os.remove(f)

def mkdir_p(dirpath):
    if not op.isdir(dirpath):
        os.mkdir(dirpath)

import sys
import pdb

class ForkedPdb(pdb.Pdb):
    """A Pdb subclass that may be used
    from a forked multiprocessing child

    """
    def interaction(self, *args, **kwargs):
        _stdin = sys.stdin
        try:
            sys.stdin = open('/dev/stdin')
            pdb.Pdb.interaction(self, *args, **kwargs)
        finally:
            sys.stdin = _stdin

##############################
#                            #
#     Homog logic            #
#                            #
##############################


def load_homog(homog_path):
    validate_file(homog_path)
    names = coord_cols + list('UXM')
    ddtype = {x: int for x in names}
    ddtype['chr'] = str
    try:
        df = pd.read_csv(homog_path, sep='\t', dtype=ddtype, names=names, comment='#')
    except:
        eprint(f'Failed loading memoization file {homog_path} . is it gzipped?')
        # os.remove(homog_path)
        return pd.DataFrame()
    return df

def pat2memfile(pat, tmp_dir_l):
    return op.join(tmp_dir_l, 'mempat_{}.homog.gz'.format(hashlib.md5(pat.encode()).hexdigest()))

def pat2homog_mp_wrap(markers, pats, tmp_dir_l, rlen, verb, force):
    # multiprocess wrapper for the pat2homog method
    # return a dict {pat: homog table}
    # homog table is markers x [U, X, M]
    params = [(markers, p, tmp_dir_l, rlen, verb, force)
               for p in pats]
    p = Pool(15)
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

def htable2str(homog_table):
    return homog_table.astype(str).agg(':'.join, axis=1)


def gen_homogs(markers, pats, tmp_dir, verb, rlen, force):

    # arguments cleanup
    markers = markers[coord_cols]          # keep only coords columns
    apats = [op.abspath(p) for p in pats]
    tmp = [validate_file(p) for p in apats]
    pats = sorted((set(apats)))            # remove duplicates, if exist
    mkdir_p(tmp_dir) # create tmp_dir if needed
    tmp_dir_l = op.join(tmp_dir, f'l{rlen}')
    mkdir_p(tmp_dir_l)

    # run pat2homog on missing markers and pats
    uxm_dict = pat2homog_mp_wrap(markers, pats, tmp_dir_l, rlen, verb, force)

    return uxm_dict


def extend_markers(markers, uniq_name):
    fake_df = markers[['chr', 'startCpG', 'endCpG']].drop_duplicates()
    fake_df['startCpG'] = (fake_df['startCpG'].astype(int) - 80).clip(lower=1)
    fake_df = fake_df.sort_values(by='startCpG')
    ext_mpath = uniq_name + '.fake_tmp_df.bed'
    fake_df.to_csv(ext_mpath, sep='\t', header=None, index=None)
    return ext_mpath

def wrap_cpp_tool(pat, markers, tmp_dir_l, rlen, verb):
    # dump extended markers file (for tabix -R)
    validate_file(pat)
    uniq_name = f'{pat2name(pat)}.tmp.{uuid.uuid4()}'
    uniq_name = op.join(tmp_dir_l, uniq_name)
    ext_mpath = extend_markers(markers, uniq_name)

    # dump reducerd markers file
    tmp_mpath = uniq_name + '.markers.bed'
    markers[coord_cols].to_csv(tmp_mpath, sep='\t', header=None, index=None)

    # homog file path
    tmp_homog = uniq_name + '.homog.gz'

    # pat to homog.gz:
    cmd = f'tabix {pat} -R {ext_mpath} | sort -k2,2n | uniq |'
    cmd += f' /cs/cbio/netanel/tools/reads_homog/homog - {uniq_name}'

    th = round(1 - (rlen - 1) / rlen, 3) + 0.001
    th1 = round((rlen - 1) / rlen, 3)
    cmd += f' -l {rlen} -r 0.0,{th},{th1},1.0 -b {tmp_mpath}'
    # cmd = f'wgbstools homog --rlen {rlen} -b {tmp_mpath} {pat} --bed'  # too slow

    so = None if verb else subprocess.PIPE
    subprocess.check_call(cmd, shell=True, stderr=so, stdout=so)
    # os.rename(op.basename(pat)[:-7] + '.uxm.bed.gz', tmp_homog)
    remove_files([ext_mpath, tmp_mpath])
    return tmp_homog


def pat2homog(markers, pat, tmp_dir_l, rlen, verb, force):

    mempat = pat2memfile(pat, tmp_dir_l)
    name = pat2name(pat)

    # load current markers if exist:

    omark = markers.copy()
    markers = markers.drop_duplicates(subset=coord_cols)
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
        eprint(msg)
        remain_mrk = markers
        homog = pd.DataFrame()

    # else, find markers not already present in mempat
    else:
        # ForkedPdb().set_trace()
        homog = load_homog(mempat)
        if homog.empty:
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
        uxm = res[['U', 'X', 'M']]
        return {pat: uxm}

    tmp_homog = wrap_cpp_tool(pat, remain_mrk, tmp_dir_l, rlen, verb)

    # homog.gz to numpy array uxm:
    uxm = load_homog(tmp_homog)
    # cleanup 
    remove_files([tmp_homog, tmp_homog + '.tbi'])
    nodump = False
    if uxm[['U', 'X', 'M']].values.sum() == 0:
        eprint('\033[91m' + f'WARNING:' + '\033[0m' +
                f' possibly failed in {pat} - all {uxm.shape[0]} ' \
                'values are zero. memoization is not updated, to be safe')
        nodump = True

    # ForkedPdb().set_trace()
    all_markers = pd.concat([homog, uxm])
    all_markers.sort_values(by=['startCpG', 'endCpG'], inplace=True)
    res = omark.merge(all_markers, how='left')
    if not nodump:
        # dump
        with gzip.open(mempat, 'wt') as mp:
            mp.write(f'# {pat}\n')
        all_markers.to_csv(mempat, sep='\t', index=None, na_rep='NA', header=None, mode='a')
        # chmod to -rw-rw-r-- so all grail users can share the same temp_dir
        os.chmod(mempat, 0o664)

    uxm = res[['U', 'X', 'M']]
    return {pat: uxm}


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
    parser.add_argument('--fast', action='store_true')
    # parser.add_argument('--threads', '-@', type=int, default=32,
            # help='number of threads to use [32]')
    parser.add_argument('--force', '-f', action='store_true',
            help='Force run homog.cpp: recreate the whole memoization table ')
    parser.add_argument('--rlen', '-l', type=int, default=-1,
            help='minimal CpGs per read required to consider the read.' \
                    ' By default this value is deduced from the atlas name (e.g. atlas.l4.tsv)')

