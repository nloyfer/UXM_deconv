#!/usr/bin/python3 -u

import os
import sys
import pandas as pd
import getpass
import os.path as op
import numpy as np
import argparse
from scipy import optimize
from homog_mem import *

def validate_ref_tissues(df, tissue_list):
    for col in tissue_list:
        if col not in df.columns:
            eprint('Invalid cell type (not in atlas):', col)
            exit()

def load_atlas(atlas_path, ignore=None, include=None):
    if not op.isfile(atlas_path):
        eprint('Invalid reference atlas (--atlas flag)')
    validate_file(atlas_path)

    # take a peek:
    df = pd.read_csv(atlas_path, sep='\t', nrows=2)
    if df.shape[1] < 8:
        eprint(f'Invalid atlas: {atlas_path}')
        exit(1)
    df = pd.read_csv(atlas_path, sep='\t')
    if df['name'].str.startswith('chr').sum() != df.shape[0]:
        eprint(f'Invalid atlas: {atlas_path}. "name" column must all start with "chr"')
        exit(1)

    if ignore is not None:
        validate_ref_tissues(df, ignore)
        for col in ignore:
            del df[col]
            df = df[df.target != col]
    elif include is not None:
        validate_ref_tissues(df, include)
        df = df[df.target.isin(include)]
        keep = list(df.columns[:8]) + include
        df = df[keep]

    df.reset_index(inplace=True, drop=True)
    return df, list(df.columns[8:])


def decon_single_samp(samp, atlas, counts, verbose, debug=False):
    """
    Deconvolve a single sample, using NNLS, to get the mixture coefficients.
    :param samp: a vector of a single sample
    :param atlas: the atlas DataFrame
    :return: the mixture coefficients
    """

    name = samp.columns[2]
    counts.columns = ['name', 'direction', 'counts']

    # remove missing sites from both sample and atlas:
    # TODO: imputation for the atlas?
    nd_cols = ['name', 'direction']
    data = samp.merge(atlas.drop_duplicates(nd_cols, ignore_index=True), on=nd_cols, how='inner').copy().dropna(axis=0)
    data = data.merge(counts.drop_duplicates(nd_cols, ignore_index=True), on=nd_cols, how='left')

    if data.empty:
        eprint(f'Warning: skipping an empty sample {name}')
        return np.nan, np.nan

    if data.shape[0] > atlas.shape[0]:
        eprint('ERROR: merge went wrong. Validate your atlas')
        return None, None
    if verbose:
        eprint('{}: {} \ {} markers'.format(name, data.shape[0], atlas.shape[0]))
    del data['name'], data['direction']

    samp = data.iloc[:, 0]
    counts = data.iloc[:, -1]
    red_atlas = data.iloc[:, 1:-1]

    # apply weights:
    red_atlas = red_atlas * counts.values[:, np.newaxis]
    samp = samp * counts

    # get the mixture coefficients by deconvolution 
    # (non-negative least squares)
    mixture, residual = optimize.nnls(red_atlas, samp)
    mixture /= np.sum(mixture)
    return mixture


def subsample(uxm, ss_rate):
    if not 1.0 >= ss_rate > 0:
        eprint('Invalid subsampling rate:', ss_rate)
        exit(1)
    if ss_rate < 1:
        uxm = np.random.binomial(uxm, ss_rate)
    return uxm

def load_pats_homog(atlas, pats, tmp_dir, verb, rlen, force, ss_rate, nodump, debug, threads):
    for pat in pats:
        if not pat.endswith('.pat.gz'):
            eprint(f'Invalid input file: {pat}. must end with .pat.gz')
            exit(1)
    pats = [op.abspath(p) for p in pats]
    uxm_dict = gen_homogs(atlas, pats, tmp_dir, verb, rlen, force, nodump, debug, threads)
    samples_df = atlas[['name', 'direction']].copy()
    counts = atlas[['name', 'direction']].copy()
    for pat in pats:
        name = pat2name(pat)
        uxm = uxm_dict[pat].values
        if ss_rate:
            uxm = subsample(uxm, ss_rate)
        counts[name] = uxm.sum(axis = 1)
        with np.errstate(divide='ignore', invalid='ignore'):
            uxm = uxm / uxm.sum(axis = 1)[:, np.newaxis]
        samples_df[name] = np.where(samples_df['direction'] == 'U', uxm[:, 0], uxm[:, 2])
    return samples_df, samples_df.columns[2:], counts


def main():
    args = parse_args()

    # load atlas:
    atlas, ref_cells = load_atlas(args.atlas, args.ignore, args.include)

    # load samples homog tables:
    rlen = deduce_l_from_name(args)
    sf, sample_names, counts = load_pats_homog(atlas, args.pats, args.tmp_dir,
                        args.verbose, rlen, args.force, args.sub_sample,
                        args.nodump, args.debug, args.threads)

    # deconvolve samples:
    df = pd.DataFrame(columns=sample_names, index=ref_cells)

    params = [(sf[['name', 'direction', samp]],
               atlas[['name', 'direction'] + ref_cells],
               counts[['name', 'direction', samp]],
               args.verbose,
               args.debug) for samp in sample_names]
    if args.debug:
        arr = [decon_single_samp(*p) for p in params]
    else:
        p = Pool(args.threads)
        arr = p.starmap(decon_single_samp, params)
        p.close()
        p.join()
    for i, samp in enumerate(sample_names):
        df[samp] = arr[i]

    # Dump results
    df.reset_index(inplace=True)
    df.columns = ['CellType'] + list(df.columns)[1:]
    df.to_csv(args.output, float_format='%.7f', index=None)
    if sys.stdout != args.output:
        eprint('dumped atlas to', args.output)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--atlas', '-a', default=DEF_ATLAS,
            help=f'path to atlas generated by build_atlas. Default: {DEF_ATLAS}')

    ref_tissues = parser.add_mutually_exclusive_group(required=False)
    ref_tissues.add_argument('--include', nargs='+',
            help='Columns in the atlas to include. Complementary to --ignore')
    ref_tissues.add_argument('--ignore', nargs='+',
            help='Columns in the atlas to remove, along with their corresponding marker lines')

    parser.add_argument('pats', nargs='+',
            help='One or more pat files to deconvolve')
    parser.add_argument('--output', '-o', default=sys.stdout,
            help='output path: a csv file. default is stdout.')
    parser.add_argument('--sub_sample', '-S', type=float,
            help='subsample from the test sample reads')
    add_memoiz_args(parser)
    return parser.parse_args()


if __name__ == '__main__':
    main()

