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
import matplotlib.pyplot as plt

print("Start procrssing.", file=sys.stderr, flush=True)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-f', '--file_dir', dest='file_dir', required=True,
        type=str,
        help='specify the file directory of pySpade process output, the Trans_genome_seq.npy file is required at this step.'
    )

    parser.add_argument(
        '-l', '--local_df', dest='local_df', required=True,
        type=str,
        help='Specify the unfiltered_local_df.csv (the output from pySpade local).'
    )

    parser.add_argument(
        '-cs', '--cutoff', dest='cutoff', required=False,
        type=float, default=0.1,
        help='Specify FDR cutoff.'
    )

    parser.add_argument(
        '-cf', '--cutoff_fc', dest='cutoff_fc', required=False,
        type=float, default=0.2, 
        help='Specify fold change cutoff, default is 20% (0.2).'
    )
    
    parser.add_argument(
        '-cx', '--cutoff_expression', dest='cutoff_expression', required=False,
        type=float, default=0.05,
        help='specify the cutoff of expressed genes. Default is 0.05 (genes expressed in more than 5 percent of cells)'
    )

    parser.add_argument(
        '-o', '--outdir', dest='outdir', required=True,
        type=str,
        help='Specify the output file directory.'
    )

    args = parser.parse_args()
    
    FILE_DIR = args.file_dir
    LOCAL_DF = args.local_df
    FDR_CUTOFF = args.cutoff
    OUTDIR = args.outdir

    #parameters
    expression_cutoff = args.cutoff_expression
    fold_change_cutoff = args.cutoff_fc
    empirical_pval_cutoff = 0.001
    hypergeom_pval_cutoff = -3
    gene_seq = np.load(FILE_DIR + 'Trans_genome_seq.npy', allow_pickle=True)
    express_level = np.load(FILE_DIR + 'Perc_cell_expr.npy')
    express_idx = np.where(express_level > expression_cutoff)[0]
    
    #read file 
    local_df = pd.read_csv(LOCAL_DF)

    express_df = local_df[local_df['gene_names'].isin(gene_seq[express_idx])\
                        &((local_df['fc_by_rand_dist_cpm'] > (1+fold_change_cutoff)) | (local_df['fc_by_rand_dist_cpm'] < (1-fold_change_cutoff)))\
                        &(local_df['Significance_score'] < FDR_CUTOFF)\
                        &(local_df['pval-empirical'] < empirical_pval_cutoff)
                        &(local_df['log(pval)-hypergeom'] < hypergeom_pval_cutoff)]


    #save the filter local df 
    express_df.to_csv(OUTDIR + '/filtered_local_df.csv')

if __name__ == '__main__':
    main()