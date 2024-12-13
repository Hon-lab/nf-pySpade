#!/usr/bin/env python3
import os
import re
import sys
import collections
import argparse
import tables
import glob
import itertools
import math

import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse
import random

print("start analyzing.", file=sys.stderr, flush=True)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-s', '--sgRNA_df', dest='sgRNA_df', required=True,
        type=str,
        help='specify the h5 file of processed sgRNA df (output from pySpade process).'
    )
    parser.add_argument(
        '-o', '--output_df', dest='output_df', required=True,
        type=str,
        help='specify the output sgrna matrix, including path and file name.'
    )
    
    args = parser.parse_args()
    sgrna_df_file = args.sgRNA_df    
    OUTPUT_DF = args.output_df
    
    sgrna_df = pd.read_hdf(sgrna_df_file, 'df')
    
    #randomized the columns and indexs
    columns = list(sgrna_df.columns)
    random.shuffle(columns)
    rows = list(sgrna_df.index)
    random.shuffle(rows)

    df_shuffled = sgrna_df.rename(columns=dict(zip(sgrna_df.columns, columns)), index=dict(zip(sgrna_df.index, rows)))

    df_shuffled.to_hdf(OUTPUT_DF, key='df', mode='w')

if __name__ == '__main__':
    main()
