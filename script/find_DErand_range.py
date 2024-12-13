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

def find_cell_num(data_dir, region):
    fc_files = glob.glob(data_dir + region + '-' + '*-foldchange')
    if len(fc_files) == 0:
        print('Missing region: ' + region, file=sys.stderr, flush=True)
        num_sgrna_cell = 0
    elif len(fc_files) == 1:
        fc_file = fc_files[0]
        num_sgrna_cell = int(fc_file.split('/')[-1].split('-')[-2])
    else:
        numbers = np.array([int(float(i.split('/')[-1].split('-')[-2])) for i in fc_files])
        chosen_num = np.max(numbers)
        fc_file = data_dir + region + '-' + str(chosen_num) + '-foldchange'
        num_sgrna_cell = int(fc_file.split('/')[-1].split('-')[-2])

    return num_sgrna_cell

def read_sgrna_dict(SGRNA_FILE):
    sgrna_dict  = {}
    with open(SGRNA_FILE) as f:
        for line in f:
            region_id, sgrna_string = line.strip().split("\t")
            sgrnas = sgrna_string.split(";")
            sgrna_dict.update({region_id : sgrnas})
    return(sgrna_dict)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-d', '--DEobs_dir', dest='DEobs_dir', required=True,
        type=str,
        help='Specify the DEobs directory, which contains the output files.'
    )

    parser.add_argument(
        '-s', '--sgrna_dict', dest='sgrna_dict', required=True,
        type=str,
        help='Specify the sgrna dictionary file.'
    )

    parser.add_argument(
        '-o', '--outdir', dest='outdir', required=True,
        type=str,
        help='Specify the output file directory.'
    )

    args = parser.parse_args()

    sgrna_dict_file = args.sgrna_dict
    DEobs_DIR = args.DEobs_dir
    OUTDIR = args.outdir

    #Read sgrna dict
    sgrna_dict = read_sgrna_dict(sgrna_dict_file)
    
    #Find the cell number of each perturbation region
    Cell_num = []
    for region in sgrna_dict.keys():
        num = find_cell_num(DEobs_DIR, region)
        Cell_num.append(num)

    #Visualize the distribution
    fig, ax = plt.subplots(figsize= (12,8))
    binwidth = 50
    ax.hist(Cell_num, 
        color='#348ABD',
        bins=range(min(Cell_num), max(Cell_num) + binwidth, binwidth))
    
    ax.set_ylabel('Frequency', fontsize=18)
    ax.set_xlabel('Number of cell', fontsize=18)
    plt.savefig(OUTDIR + '/Cell_num_distribution.pdf')

    #save to output
    if len(Cell_num) < 10:
        output_file = open(OUTDIR + '/bin.txt', 'w')
        for i in Cell_num:
            output_file.write(str(i) + '\n')

    else:
        output_file = open(OUTDIR + '/bin.txt', 'w')
        output_file.write(str(round(np.quantile(Cell_num, 0.1))) + '\n')
        output_file.write(str(round(np.quantile(Cell_num, 0.25))) + '\n')
        output_file.write(str(round(np.quantile(Cell_num, 0.4))) + '\n')
        output_file.write(str(round(np.quantile(Cell_num, 0.5))) + '\n')
        output_file.write(str(round(np.quantile(Cell_num, 0.6))) + '\n')
        output_file.write(str(round(np.quantile(Cell_num, 0.75))) + '\n')
        output_file.write(str(round(np.quantile(Cell_num, 0.9))) + '\n')

    output_file.close()

if __name__ == '__main__':
    main()