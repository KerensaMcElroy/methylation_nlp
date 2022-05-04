#importing libraries
from tkinter import W
import torch
import numpy as np
import random
import matplotlib.pyplot as plt
import sklearn.metrics as metrics

torch.manual_seed(1)
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f'Using {device} device')

#defining dataset class
from torch.utils.data import Dataset, DataLoader, random_split
import pandas as pd
from collections import Counter

class DnaMethylation(torch.utils.data.Dataset):
    def __init__(self, data_file):
        self.df = pd.read_csv(data_file,sep='\t',header=None)
        self.oligos = self.df[0]
        self.uniq_oligos = self.get_uniq_oligos()
        self.oligo_to_index = {oligo: index for index, oligo in enumerate(self.uniq_oligos)}
        self.meth = np.array(self.df[1].str.split(',', expand=True).astype(float))

    def get_uniq_oligos(self):
        all_oligos = self.oligos.str.cat(sep=',').split(',')
        oligo_counts = Counter(all_oligos)
        return sorted(oligo_counts, key=oligo_counts.get, reverse=True)
    
    def balance(self, binary_meth):
        weights = np.copy(binary_meth)
        num_pos = sum(binary_meth)
        zero_indices = np.where(binary_meth == 0)[0].tolist()
        zero_subset = random.sample(zero_indices, num_pos)
        weights[zero_subset] = 1
        return weights
    
    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):

        oligo_list = self.df[0][idx].split(',')
        oligo_indexes = [self.oligo_to_index[o] for o in oligo_list]
        binary_meth = np.where(self.meth[idx] > 0, 1, 0)
        weight = self.balance(binary_meth)
        return (
            torch.tensor(oligo_indexes),
            torch.tensor(binary_meth),
            torch.tensor(weight)
        )

#parameters
b_size = 64

#create dataloader
meth_file = 'window_mean_HE1_chrm11_iaN7_win10_step10_frag1000.txt'
dataset = DnaMethylation(meth_file)
train_loader = DataLoader(dataset, batch_size=b_size, shuffle=True)

#defining the network
from torch import nn
import torch.nn.functional as F

class LSTMTagger(nn.Module):

    def __init__(self, embedding_dim, hidden_dim, relu_dim, oligo_set_size):
        super(LSTMTagger, self).__init__()
        self.hidden_dim = hidden_dim
        self.oligo_embeddings = nn.Embedding(oligo_set_size, embedding_dim)

        # The LSTM takes oligonucleotide embeddings as inputs, and outputs hidden states
        # with dimensionality hidden_dim.
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, batch_first=True)

        # The linear layer that maps from hidden state space to tag space
        #self.hidden2tag = nn.Linear(hidden_dim, 2)
        self.linear_relu_stack = nn.Sequential(     #trying a relu stack instead
            nn.Linear(hidden_dim, relu_dim),
            nn.ReLU(),
            nn.Linear(relu_dim, relu_dim),
            nn.ReLU(),
            nn.Linear(relu_dim,1),
        )
       # self.sigmoid = nn.Sigmoid()
        x = self.flatten(nuc_lstm_out[-1])
        tag_space = self.linear_relu_stack(x)
        tag_scores = F.log_softmax(tag_space, dim=1)

        return tag_scores

    def forward(self, fragment):
       # print('fragment shape:',fragment.shape)
        embeds = self.oligo_embeddings(fragment)
       # print(len(fragment))
       # print('embed shape:',embeds.shape)
        lstm_out, _ = self.lstm(embeds)
       # print('lstm shape:', lstm_out.shape)
        linear_out = self.linear_relu_stack(lstm_out)
        probs = self.sigmoid(linear_out)

       # print('tagspace shape:' ,tag_space.shape)
        #tag_scores = F.log_softmax(tag_space, dim=1)
        #print('tag_scores shape:',tag_scores.shape)
        #print('tagscores:',tag_scores)
        return probs

#hyper parameters
embedding_dim = 20
hidden_dim = 200
relu_dim = 10
learning_rate = 0.01
epochs = 1000



for j,(fragment_train,meth_train, weight) in enumerate(train_loader):
   # print(j, fragment_train, meth_train, weight)
    print(torch.sum(meth_train, dim=1))
    print(torch.sum(weight, dim=1))


