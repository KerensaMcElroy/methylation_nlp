import torch
from torch.utils.data import Dataset
from torch.utils.data import DataLoader, random_split

import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np

torch.manual_seed(1)
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f'Using {device} device')

#code to make dataset
#data in file looks like:
#
#8819 8869 0.000 CGAGAGGAACACTGAGGTTCCTGGCACCACTTCCTCTGAGCCCCATCTCC
#8844 8894 0.017 ACCACTTCCTCTGAGCCCCATCTCCCCTCCTGATCTGGACAGGAGGGTCG
#8869 8919 0.017 CCTCCTGATCTGGACAGGAGGGTCGATTCCCCTGCCTTGTCTGGAAGGGG
#8894 8944 0.028 ATTCCCCTGCCTTGTCTGGAAGGGGTTCACGACCTTCCCGTGGCACCTCA
#
#columns are start, end, mean methylation, sequence
class MethDataset(Dataset):

    def __init__(self, src_file, num_rows=None):
        self.o_data = np.loadtxt(src_file, max_rows=num_rows,
            usecols=3, delimiter=" ",
            skiprows=0, dtype='U50')
        self.m_data = np.loadtxt(src_file, max_rows=num_rows,
            usecols=2, delimiter=" ", skiprows=0,
            dtype=np.float32)
        self.s_data = np.loadtxt(src_file, max_rows=num_rows,
            usecols=0, delimiter=" ",
            skiprows=0, dtype=np.int32)
        self.e_data = np.loadtxt(src_file, max_rows=num_rows,
            usecols=0, delimiter=" ",
            skiprows=0, dtype=np.int32) 

    def __len__(self):
        return len(self.m_data)

    def __getitem__(self, idx):
        oligo = self.o_data[idx]
        meth = self.m_data[idx]
        start = self.s_data[idx]
        end = self.e_data[idx]
    
        return (oligo, meth, start, end)




if __name__ == '__main__':
    dataset = MethDataset('/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/users/mce03b/methylation_nlp/analysis/window_mean_HE1_chrm2_iaN7_win50_step25.txt')
    trainset, valset = random_split(dataset, [round(len(dataset)*0.7), round(len(dataset)*0.3)])
    train_loader = DataLoader(trainset, batch_size=10, shuffle=True, num_workers=2)
    val_loader = DataLoader(valset, batch_size=10, shuffle=True, num_workers=2)

    for i, batch in enumerate(train_loader):
        print(i, batch)

    for i, batch in enumerate(val_loader):
        print(i, batch)
