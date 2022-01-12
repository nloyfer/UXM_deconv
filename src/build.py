#!/usr/bin/python3 -u

import os
import sys
import pandas as pd
import numpy as np
import os.path as op
import argparse
import subprocess
from multiprocessing import Pool
from homog_mem import *


##############################
#                            #
#     Loading stuff          #
#                            #
##############################

def load_markers(mpath, use_um):
    validate_file(mpath)
    names = coord_cols + ['target', 'name']
    cols = list(range(7))
    if not use_um:
        names += ['mt', 'mb']
        cols += [9, 10]
    # find out if there is a header:
    peek = pd.read_csv(mpath, sep='\t', header=None, nrows=1)
    skiprows = 0 if str(peek.iloc[0, 2]).isdigit() else 1
    df = pd.read_csv(mpath, sep='\t', names=names, usecols=cols, skiprows=skiprows)

    # validation: make sure name column is unique for each marker:
    if df['name'].str.startswith('chr').sum() != df.shape[0]:  # todo: some references don't start with "chr"!
        eprint(f'Invalid marker file: {mpath}. "name" column must start with "chr"')
        exit(1)
    df['direction'] = 'U'
    if use_um:
        tf = df.copy()
        tf['direction'] = 'M'
        df = pd.concat([df, tf]).reset_index(drop=True)
    else:
        df.loc[df['mt'] > df['mb'], 'direction'] = 'M'
        del df['mt'], df['mb']
    return df


def load_groups(groups_path, pats):
    pats = [op.abspath(validate_file(p)) for p in pats]
    pats = drop_dup_keep_order(pats)
    if not groups_path:
        groups = [pat2name(p) for p in pats]
        if len(set(groups)) != len(pats):
            eprint('[uxm build] Error: some pats have identical names!')
            exit(1)  # TODO: use custom exceptions
        return pd.DataFrame({'group': groups, 'full_path': pats})

    validate_file(groups_path)
    df = pd.read_csv(groups_path, comment='#', index_col=False)
    if 'include' in df.columns:
        df = df[df['include']]
    df = df[['name', 'group']]
    # TODO: document this, or remove it (but check groups are mutually exclusive?)
    df = df[~df['group'].str.contains(':')] # remove merged groups. 
    df['full_path'] = df['name'].replace(ref_pat_full_paths(df, pats))
    return df


def ref_pat_full_paths(groups_df, pats): # TODO: will fail with name subset another name!
    nd = {}
    for name in set(groups_df['name']):

        files = [f for f in pats if name
                in f and f.endswith('.pat.gz')]
        if ':' not in name:
            files = [f for f in files if ':' not in f]  # TODO: address the ":" issue
        if not files:
            eprint(f'Error: no {name} found in supplied pat files')
            exit(1)
        if len(files) > 1:
            eprint(f'Error: ambiguous name {name} in pats')
            eprint('\n'.join(files))
            exit(1)
        res = files[0]
        nd[name] = res
    return nd


##############################
#                            #
#     Main                   #
#                            #
##############################


def add_col(df, uxm_dict, groups_df, group):
    uxm = np.zeros((df.shape[0], 3), dtype=float)
    for m in groups_df[groups_df['group'] == group]['full_path']:
        uxm += uxm_dict[m].values
    with np.errstate(divide='ignore', invalid='ignore'):
        uxm = uxm / uxm.sum(axis = 1)[:, np.newaxis]
    res = np.where(df['direction'] == 'U', uxm[:, 0], uxm[:, 2])
    df[group] = np.round(res, 3)


def main():
    args = parse_args()
    # load markers and groups files
    df = load_markers(args.markers, args.use_um)
    groups_df = load_groups(args.groups, args.pats)

    # calc homog tables:
    uxm_dict = gen_homogs(df, groups_df['full_path'], args.tmp_dir,
                          args.verbose, args.rlen, args.force, args.threads)

    # merge to groups and dump
    for group in sorted(groups_df['group'].unique()):
        add_col(df, uxm_dict, groups_df, group)

    # dump atlas
    df.to_csv(args.output, sep='\t', index=None, na_rep='NA')
    eprint('dumped atlas:', args.output)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--use_um', action='store_true',
            help='Use U and M values for all markers. Otherwise deduce if the marker ' \
            'is U or M from the markers file')
    parser.add_argument('--markers', '-m', required=True,
            help='path to markers file')
    parser.add_argument('--groups', '-g',
            help=f'path to group file. If none given, '\
                   ' each pat gets its own column in the reference atlas')
    parser.add_argument('--output', '-o', required=True,
            help='output path for the atlas')
    parser.add_argument('--pats', required=True, nargs='+',
            help=f'Reference pat files.') # todo: 1. beta hack. 2. Document it somewhere (groups, pats, subset)
    add_memoiz_args(parser)
    return parser.parse_args()


if __name__ == '__main__':
    main()

