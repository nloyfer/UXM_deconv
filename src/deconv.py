#!/usr/bin/python3 -u

import os
import sys
import pandas as pd
import getpass
import os.path as op
import numpy as np
import argparse
from scipy import optimize
# from plot_dec_res import plot_res
from homog_mem import *


def load_atlas(atlas_path, ignore=None):
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
        for col in ignore:
            if col not in df.columns:
                eprint('Invalid cell type (not in atlas):', col)
                exit()
            del df[col]
            df = df[df.target != col]
    df.reset_index(inplace=True, drop=True)
    return df, list(df.columns[8:])


def decon_single_samp(samp, atlas, counts, weighted, verbose, debug=False):
    """
    Deconvolve a single sample, using NNLS, to get the mixture coefficients.
    :param samp: a vector of a single sample
    :param atlas: the atlas DataFrame
    :return: the mixture coefficients
    """

    name = samp.columns[1]
    counts.columns = ['name', 'direction', 'counts']

    # remove missing sites from both sample and atlas:
    # TODO: imputation for the atlas?
    data = samp.merge(atlas.drop_duplicates(['name', 'direction'], ignore_index=True), on=['name', 'direction'], how='inner').copy().dropna(axis=0)
    data = data.merge(counts.drop_duplicates(['name', 'direction'], ignore_index=True), on=['name', 'direction'], how='left')

    if data.empty:
        eprint(f'Warning: skipping an empty sample {name}')
        return np.nan, np.nan

    if data.shape[0] > atlas.shape[0]:
        eprint('ERROR: merge went wrong. Validate your atlas')
        return None, None
    if verbose:
        eprint('{}: {} \ {} sites'.format(name, data.shape[0], atlas.shape[0]))
    del data['name'], data['direction']

    samp = data.iloc[:, 0]
    counts = data.iloc[:, -1]
    red_atlas = data.iloc[:, 1:-1]
    if weighted:
        red_atlas = red_atlas * counts.values[:, np.newaxis]
        samp = samp * counts

    # get the mixture coefficients by deconvolution 
    # (non-negative least squares)
    mixture, residual = optimize.nnls(red_atlas, samp)
    mixture /= np.sum(mixture)
    return mixture, residual


def subsample(uxm, ss_rate):
    if not 1.0 >= ss_rate > 0:
        eprint('Invalid subsampling rate:', ss_rate)
        exit(1)
    if ss_rate < 1:
        uxm = np.random.binomial(uxm, ss_rate)
    return uxm

def load_pats_homog(atlas, pats, tmp_dir, verb, rlen, force, ss_rate, threads):
    pats = [op.abspath(p) for p in pats]
    uxm_dict = gen_homogs(atlas, pats, tmp_dir, verb, rlen, force, threads)
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


def dump_counts(args, counts, atlas):
    if getpass.getuser() != 'nloyfer':
        return
    counts.loc[atlas.iloc[:,3:].isnull().any(axis=1).values] = np.nan
    counts = pd.concat([atlas[['name', 'target']], counts], axis=1)
    counts.to_csv(args.prefix + '.count.csv', index=None, float_format='%.0f')


def main():
    args = parse_args()

    # load atlas:
    atlas, ref_cells = load_atlas(args.atlas, args.ignore)

    # load samples homog tables:
    rlen = deduce_l_from_name(args)
    sf, sample_names, counts = load_pats_homog(atlas, args.pats, args.tmp_dir,
                        args.verbose, rlen, args.force, args.sub_sample, args.threads)

    # deconvolve samples:
    df = pd.DataFrame(columns=sample_names, index=ref_cells)
    rf = pd.DataFrame(columns=['Residuals'], index=sample_names)

    params = [(sf[['name', 'direction', samp]],
               atlas[['name', 'direction'] + ref_cells],
               counts[['name', 'direction', samp]],
               args.weighted,
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
        df[samp], rf.loc[samp] = arr[i]

    # Dump results
    dec_out = args.prefix + '.dec_out.csv'
    df.to_csv(dec_out, float_format='%.7f')
    if getpass.getuser() == 'nloyfer':
        f = dec_out
        t = f + '.tmp'
        os.system(f'cat {f} | change_float_format -s , -d 6 > {t} && mv -f {t} {f}')
    if args.resid:
        rf.to_csv(args.prefix + '.residuals.csv', float_format='%.6f')
    if args.counts:
        dump_counts(args, counts.iloc[:, 1:], atlas[['name', 'target', 'direction'] + ref_cells])

    # Plot pie charts
    # if args.plot or args.png:
        # plot_res(df, args.prefix, args.plot, args.plot_stubs,
                 # not args.plot_nocollapse, args.plot_full_legend)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--atlas', '-a', required=True,
            help='path to atlas generated by build_atlas')
    parser.add_argument('--ignore', nargs='+',
            help='Columns in the atlas to remove, along with their corresponding marker lines')
    parser.add_argument('pats', nargs='+',
            help='One or more pat files to deconvolve')
    parser.add_argument('--prefix', '-p', default='./out',
            help='prefix for output files (csv and png)')
    parser.add_argument('--resid', '-r', action='store_true',
            help='If set, output residuals of the nnls')
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('--weighted', '-w', action='store_true',
            help='weight each marker by the number of reads it contains in the input sample (recommended)')
    parser.add_argument('--counts', action='store_true',
            help='output read-count table [nr_blocks x nr_input_samples]')
    parser.add_argument('--sub_sample', '-S', type=float,
            help='subsample from the test sample reads')
    add_memoiz_args(parser)
    return parser.parse_args()


if __name__ == '__main__':
    main()

