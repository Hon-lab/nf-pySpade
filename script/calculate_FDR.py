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
        '-g', '--global_df', dest='global_df', required=True,
        type=str,
        help='Specify the unfiltered_global_df.csv (the output from pySpade global).'
    )

    parser.add_argument(
        '-r', '--randomized_global_df', dest='randomized_global_df', required=True,
        type=str,
        help='Specify the unfiltered_global_df.csv from randomized sgRNA matrix.'
    )

    parser.add_argument(
        '-c', '--cutoff', dest='cutoff', required=False,
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
    GLOBAL_DF = args.global_df
    FP_GLOBAL = args.randomized_global_df
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
    global_df = pd.read_csv(GLOBAL_DF)

    cutoff_list = list(-np.arange(1, 50, 0.5))
    num_hits = []
    for Significance_score_cutoff in cutoff_list:
        express_df = global_df[global_df['gene_names'].isin(gene_seq[express_idx])\
                            &((global_df['fc_by_rand_dist_cpm'] > (1+fold_change_cutoff)) | (global_df['fc_by_rand_dist_cpm'] < (1-fold_change_cutoff)))\
                            &(global_df['Significance_score'] < Significance_score_cutoff)\
                            &(global_df['pval-empirical'] < empirical_pval_cutoff)\
                            &(global_df['log(pval)-hypergeom'] < hypergeom_pval_cutoff)]
        num, _ = express_df.shape
        num_hits.append(num)


    #false positive
    fp_global = pd.read_csv(FP_GLOBAL)
    num_false_hits = []
    for Significance_score_cutoff in cutoff_list:
        express_df = fp_global[fp_global['gene_names'].isin(gene_seq[express_idx])\
                            &((fp_global['fc_by_rand_dist_cpm'] > (1+fold_change_cutoff)) | (fp_global['fc_by_rand_dist_cpm'] < (1-fold_change_cutoff)))\
                            &(fp_global['Significance_score'] < Significance_score_cutoff)\
                            &(fp_global['pval-empirical'] < empirical_pval_cutoff)\
                            &(fp_global['log(pval)-hypergeom'] < hypergeom_pval_cutoff)]
        num, _ = express_df.shape
        num_false_hits.append(num)

    #calculate FDR
    false_pos_both = np.array(num_false_hits)
    true_both = np.array(num_hits)
    FDR = np.array(false_pos_both / (true_both + false_pos_both))

    #plot FDR
    fig, ax = plt.subplots(figsize=(8,6))
    ax.grid(ls='--', color='#D8D8D8')
    ax.set_xlabel('-Significance score cutoff', fontsize=18)
    ax.set_ylabel('FDR', fontsize=18)
    ax.set_xticks(np.arange(0, 50, 1))
    ax.set_xticklabels(np.arange(0,50,1), rotation=-45)
    ax.set_yticks(np.arange(0, 1, 0.05))
    ax = plt.plot(np.array(np.arange(1, 50, 0.5)),
                FDR, 
                color='#348ABD', marker='o', markersize=6)
    plt.savefig(OUTDIR + '/FDR.pdf')

    #get the cutoff 
    if np.sum(FDR < FDR_CUTOFF) == 0: 
        print("No Significance score cutoff meets the FDR threshold. Please refer to FDR.pdf file to view the FDR curve.", file=sys.stderr, flush=True)
    else:
        Significance_cutoff = cutoff_list[np.argwhere(FDR < FDR_CUTOFF)[0][0]]
        output_file = open(OUTDIR + '/Significance_score_cutoff_FDR.txt', 'w')
        output_file.write(str(Significance_cutoff) + '\n')
        output_file.close()

if __name__ == '__main__':
    main()