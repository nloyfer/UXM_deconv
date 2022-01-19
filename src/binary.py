#!/usr/bin/python3 -u

import os
import sys
import pandas as pd
import getpass
import os.path as op
import numpy as np
import argparse
from homog_mem import *
from build import load_markers
from deconv import load_pats_homog



def decon_single_samp(samp, counts, weighted, verbose, debug=False):
    """
    Deconvolve a single sample, using NNLS, to get the mixture coefficients.
    :param samp: a vector of a single sample
    :param atlas: the atlas DataFrame
    :return: the mixture coefficients
    """

    name = samp.columns[2]
    counts.columns = ['name', 'direction', 'counts']

    # remove missing sites from both sample and atlas:
    nd_cols = ['name', 'direction']
    data = samp.merge(counts.drop_duplicates(nd_cols, ignore_index=True), on=nd_cols, how='left').copy().dropna(axis=0)

    if data.empty:
        eprint(f'Warning: skipping an empty sample {name}')
        return np.nan, np.nan

    if verbose:
        eprint('{}: {} \ {} sites'.format(name, data.shape[0], samp.shape[0]))
    del data['name'], data['direction']

    samp = data.iloc[:, 0]
    counts = data.iloc[:, -1]
    if not weighted:
        return samp.mean()
    return np.average(samp, weights=counts)


def load_binary_markers(mpath):
    markers = load_markers(mpath, False)
    targets = markers['target'].unique().tolist()
    if len(targets) > 1:
        eprint(f'[uxm]: WARNING: more than one targets: {targets}')
    target = targets[0]
    return markers, target

def main():
    args = parse_args()

    # load markers and target:
    markers, target = load_binary_markers(args.markers)

    rlen = args.rlen
    if rlen < 0:
        rlen = 3
    # load samples homog tables:
    sf, sample_names, counts = load_pats_homog(markers, args.pats, args.tmp_dir,
                        args.verbose, rlen, args.force, args.sub_sample, args.threads)

    # deconvolve samples:
    df = pd.DataFrame(columns=sample_names, index=[target])

    params = [(sf[['name', 'direction', samp]],
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
        df[samp] = arr[i]

    # Dump results
    dec_out = args.prefix + '.dec_out.csv'
    if args.transpose:
        df = df.T
    df.to_csv(dec_out, float_format='%.7f')
    if getpass.getuser() == 'nloyfer':
        f = dec_out
        t = f + '.tmp'
        os.system(f'cat {f} | change_float_format -s , -d 6 > {t} && mv -f {t} {f}')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--markers', '-m', required=True,
            help='path to markers file')
    parser.add_argument('pats', nargs='+',
            help='One or more pat files to deconvolve')
    parser.add_argument('--prefix', '-p', default='./out',
            help='prefix for output files (csv and png)')
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('--weighted', '-w', action='store_true',
            help='weight each marker by the number of reads it contains in the input sample (recommended)')
    parser.add_argument('--sub_sample', '-S', type=float,
            help='subsample from the test sample reads')
    parser.add_argument('--transpose', action='store_true',
            help='transpose output')
    add_memoiz_args(parser)
    return parser.parse_args()


if __name__ == '__main__':
    main()


