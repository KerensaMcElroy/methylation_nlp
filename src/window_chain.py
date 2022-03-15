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

#global logging
logging.basicConfig(filename = 'window_func.log', 
                    level=logging.INFO, 
                    filemode='w', 
                    format='%(asctime)s - %(levelname)s - %(message)s')

logger = logging.getLogger()


def make_frag(in_file, frag_len):

    oligos=[]
    methyl=[]
    frag=0
    
    stem, ext = in_file.split('.')
    out_file = stem+'_frag'+str(frag_len)+'.'+ext

    with open(in_file) as in_f:
        with open(out_file, 'w') as out_f:
            line = in_f.readline()
            while line:
                if frag < frag_len:
                    parts = line.split()
                    seq = parts[-1]
                    meth = parts[-2]
                    oligos.append(seq)
                    methyl.append(meth)
                    frag+=len(seq)
                    line = in_f.readline()
                else:
                    data = ','.join(oligos)+'\t'+','.join(methyl)+'\n'
                    out_f.write(data)
                    oligos=[]
                    methyl=[]
                    frag = 0

def main(args):
    """Main logic of program."""
    frag_len = args.frag
    in_file = args.in_file

    make_frag(in_file, frag_len)
    
if __name__ == "__main__":

    import argparse


    #parse command line arguments
    parser = argparse.ArgumentParser(description='Converts sliding windows from window_func.py into \
        longer fragments with multiple tags, mimicking sentence / tagged word structure for lstm')
                                     
    #positional arguments
    parser.add_argument('frag', help='Length of fragment', type=int)
    parser.add_argument('in_file', help='Output file of window_func.py', type=str)

    args = parser.parse_args()
    #run program
    main(args)
