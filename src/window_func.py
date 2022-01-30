#!/usr/bin/env python

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import time
import numpy as np
import pandas as pd
import math
import os
import logging
import gzip
from Bio import SeqIO

<<<<<<< HEAD
#global logging
logging.basicConfig(filename = 'window_func.log', 
                    level=logging.INFO, 
                    filemode='w', 
                    format='%(asctime)s - %(levelname)s - %(message)s')

logger = logging.getLogger()


def make_windows(out_file, genome, win, slide, stress,chrom,animal):

    df = pd.read_csv('/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/Data_for_Tai/filtered_meth_HE%s_chrm%s.dat'%(stress,chrom), sep="\s+")
    logger.info("count data shape: %s", df.shape)

    n_names = [col for col in df.columns if col.startswith('N')]
    logger.info("n_names: %s", n_names)

    x_names = [col for col in df.columns if col.startswith('X')]
    logger.info("x_names: %s", x_names)
=======
global logger
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
fh = logging.FileHandler('window_fun.log')
fh.setLevel(logging.INFO)
log.addHandler(fh)

def fasta_slider(seq_record, size, step):
    """generates sliding windows over nucleotide sequence"""
    i = 0
    while (i < len(seq_record.seq)) :
        f_out = (str(seq_record.seq[i:i+size]) + "\n")
        i += step
    
def main(args):
    """Main logic of program."""
    genome = args.genome
    win = args.win
    slide = args.slide
    stress = args.HE
    chrom = args.chrom
    animal = args.an

    df = pd.read_csv('test/data/filtered_meth_HE%s_chrm%s.dat'%(stress,chrom), sep="\s+")
    log.info("count data shape:", df.shape)

    n_names = [col for col in df.columns if col.startswith('N')]
    log.info("n_names:", n_names)

    x_names = [col for col in df.columns if col.startswith('X')]
    log.info("x_names:", x_names)
>>>>>>> 62fda3804f45e591b953be7df4689c368d04bf02

    n_pos = df.shape[0]
    n_ani = len(n_names)

    a = np.zeros((n_pos, n_ani))
    for i in range(n_ani):
        a[:,i] = df[x_names[i]].values/df[n_names[i]].values

    ## convert NaN to 0
    a = np.nan_to_num(a)
<<<<<<< HEAD
    logger.info("a.shape: %s", a.shape)
=======
    log.info("a.shape:", a.shape)
>>>>>>> 62fda3804f45e591b953be7df4689c368d04bf02

##================================================================================================
    pos = df["Pos"]
    site_end = int(pos.max())
<<<<<<< HEAD
=======
    print("site_end:", site_end)
>>>>>>> 62fda3804f45e591b953be7df4689c368d04bf02


## sliding windows for single animal
    ani_dict={}
    for j, n in enumerate(n_names):
        ani_dict[n]=j
    i_a=ani_dict[animal]
<<<<<<< HEAD
    logger.debug('Animal id: %s, animal index %s', animal, i_a)
=======
    log.debug('Animal id: %s, animal index %s', animal, i_a)
>>>>>>> 62fda3804f45e591b953be7df4689c368d04bf02
    start_time = time.time()


#impute zeros to empty sites
    y = np.zeros(site_end) # +1 because from 0 to df[:,0].max()
    for i,j in enumerate(pos):
<<<<<<< HEAD
=======
        print(i,j,i_a)
        print(a[i,i_a])
>>>>>>> 62fda3804f45e591b953be7df4689c368d04bf02
        y[int(j)-1] = a[i,i_a]
# find where the zeros end
    site_start = int(pos.min())

# calculate window_mean
    n_windows = math.ceil((len(y[site_start:]) - win)/slide) + 1
    window_mean = np.zeros(n_windows)
    window_seq = np.empty(n_windows, dtype='S200') #max length window possible is 200
    window_beg = np.zeros(n_windows)
    window_end = np.zeros(n_windows)
    i = 0 #initiate indexing for window_mean
    p = site_start-1 #because python uses 0 indexing
    record_dict = SeqIO.index(genome,'fasta')
    record = record_dict[chrom]
    while i < n_windows:
        w_mean = y[p:p+win].mean()
<<<<<<< HEAD
=======
        print(y[p:p+win])
>>>>>>> 62fda3804f45e591b953be7df4689c368d04bf02
        window_mean[i]=w_mean
        window_beg[i]=p
        window_end[i]=window_beg[i]+win
        window_seq[i]=str(record.seq[int(window_beg[i]):int(window_end[i])])
<<<<<<< HEAD
=======
        print(window_seq[i]) 
>>>>>>> 62fda3804f45e591b953be7df4689c368d04bf02
        p += slide
        i+=1 

    #print("finished -- ia: {}, time: {:.2f}".format(animal, time.time() - start_time))
    ta=np.array(list(zip(window_beg+1,window_end+1,window_mean,window_seq)), 'd,d,f,U200')
<<<<<<< HEAD
    np.savetxt(out_file, ta, fmt="%d %d %.3f %s")

def get_out_file_name(stress, chrom, animal, win, slide):  
    return "window_mean_HE%s_chrm%s_ia%s_win%s_step%s.txt"%(stress, chrom, animal, win, slide)

#unused for now
=======
    np.savetxt("window_mean_HE%s_chrm%s_ia%s_win%s_step%s.txt"%(stress, chrom, animal, win, slide),ta , fmt="%d %d %.3f %s")

>>>>>>> 62fda3804f45e591b953be7df4689c368d04bf02
def fasta_slider(seq_record, size, step):
    """generates sliding windows over nucleotide sequence"""
    i = 0
    while (i < len(seq_record.seq)) :
        f_out = (str(seq_record.seq[i:i+size]) + "\n")
        i += step
<<<<<<< HEAD
    
def main(args):
    """Main logic of program."""
    genome = args.genome
    win = args.win
    slide = args.slide
    stress = args.HE
    chrom = args.chrom
    animal = args.an

    out_file = get_out_file_name(stress, chrom, animal, win, slide)

    make_windows(out_file, genome, win, slide, stress, chrom, animal)
    
=======

>>>>>>> 62fda3804f45e591b953be7df4689c368d04bf02
if __name__ == "__main__":

    import argparse

<<<<<<< HEAD

    #parse command line arguments
    parser = argparse.ArgumentParser(description='Generates tab-delimited tidy\
        data for sliding windows over methylation input files generated by Andrew.')
=======
# setup logging
    logging.basicConfig(filename='./window.log', level=logging.INFO,
                        format='%(asctime)s %(messages)s')

    #parse command line arguments
    parser = argparse.ArgumentParser(description='Generates tab-delimited tidy\
        data for sliding windows over methylation input files generated by xx.')
>>>>>>> 62fda3804f45e591b953be7df4689c368d04bf02
                                     
    #positional arguments
    parser.add_argument('genome', help='Input fasta genome file')
    parser.add_argument('win', help='Sliding window size', type=int)
    parser.add_argument('slide', help='Sliding window step', type=int)
    parser.add_argument("HE", help="Heat stress stage, one of [1,2,3]")
    parser.add_argument("chrom", help="Chromosome id, from 1 to 29")
<<<<<<< HEAD
    parser.add_argument("an", help="Animal id, one of [N1,N4,N7,N10,N13,N16,N19,N22,N25,N28,N31,N34,N40,Nj43]", type=str)

    args = parser.parse_args()
=======
    parser.add_argument("an", help="Animal id, one of [2,5,8,11,14,17,20,23,26,29,32,35,41,44]", type=str)
#    parser.add_argument('in_meth', help='Input methylation file')

    args = parser.parse_args()

>>>>>>> 62fda3804f45e591b953be7df4689c368d04bf02
    #run program
    main(args)
