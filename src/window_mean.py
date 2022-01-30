#!/usr/bin/env python
import numpy as np
import pandas as pd
import sys
import os
import logging
import time
#import matplotlib.pyplot as plt
##================================================================================================

#global logger
#log = logging.getLogger(__name__)


## Parameters:
np.random.seed(1)
window_size = 50

#i_chrm = 1
#ia_min = 1

## heat stress state
i_HE = int(sys.argv[1])
print("i_HE:", i_HE)

## chrm id
i_chrm = int(sys.argv[2])
print("i_chrm:", i_chrm)

## animal id
ia_min = int(sys.argv[3])
print("ia_min:", ia_min)

#ia_step = int(sys.argv[4])
#print("ia_step:", ia_step)


df = pd.read_csv('../Data_for_Tai/filtered_meth_HE%s_chrm%s.dat'%(i_HE,i_chrm), sep="\s+")
print("df.shape:", df.shape)

n_names = [col for col in df.columns if col.startswith('N')]
print("n_names:", n_names)
print("len(n_names):", len(n_names))

x_names = [col for col in df.columns if col.startswith('X')]
print("x_names:", x_names)
print("len(x_names):", len(x_names))

##================================================================================================
n_pos = df.shape[0]
n_ani = len(n_names)

a = np.zeros((n_pos, n_ani))
for i in range(n_ani):
    a[:,i] = df[x_names[i]].values/df[n_names[i]].values
    
## convert NaN to 0
a = np.nan_to_num(a)
print("a.shape:", a.shape)

##================================================================================================
pos = df["Pos"]
site_end = int(pos.max())
print("site_end:", site_end)

##--------------
#for ia in range(ia_min, ia_min + ia_step):
ia = ia_min

print("ia:", ia)
start_time = time.time()

# impute zeros to empty sites
y = np.zeros(site_end+1) # +1 because from 0 to df[:,0].max()
for i,j in enumerate(pos):
    y[int(j)] = a[i,ia]
    
# cut the beginning-zero part
site_start = int(pos.min())
y = y[site_start:]

# calculate window_mean
n_windows = len(y) - window_size + 1
window_mean = np.zeros(n_windows)
for i in range(n_windows):
    window_mean[i] = y[i:i+window_size].mean()
    #print(i, y[i:i+window_size], window_mean[i])

#bins = np.linspace(0,1,21)
#print("bins:", bins)
#hist, bin_edges = np.histogram(window_mean, bins)
#np.savetxt("hist_HE%s_chrm%s_ia%s.txt"%(i_HE, i_chrm, ia), np.array((bin_edges[:-1], hist)).T, fmt="%f %f")
#print("finished -- ia: {}, time: {:.2f}".format(ia, time.time() - start_time))

np.savetxt("window_mean_HE%s_chrm%s_ia%s.txt"%(i_HE, i_chrm, ia), window_mean, fmt="%f")
#print("finished -- ia: {}, time: {:.2f}".format(ia, time.time() - start_time))

def main(args):
    """Main logic of program."""

if __name__ == "__main__":

    import argparse

    #set up logging
    #logging.basicConfig(filename='./vcf_to_structure.log', level=logging.INFO,
    #                    format='%(asctime)s %(messages)s')

    #parse command line arguments
    parser = argparse.ArgumentParser(description='Converts a vcf file output\
                                     by GATK to a structure input format')
    #positional arguments
    parser.add_argument('in_file', help='Input vcf file')

    args = parser.parse_args()

    #run program
    main(args)
