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
        zero_indicies = np.where(binary_meth == 0)[0].tolist()
        if  len(zero_indicies) >= num_pos:
            zero_subset = random.sample(zero_indicies, num_pos)
            weights[zero_subset] = 1
        else:
            weights[zero_indicies] = 1
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

#create dataloader - only 1000bp windows with at least one meth site
meth_file = 'window_mean_HE1_chrm11_iaN7_win10_step10_frag1000.txt'
dataset = DnaMethylation(meth_file)
counts = torch.tensor(np.sum(dataset.meth, axis = 1))
one_indices = torch.nonzero(counts).flatten().tolist()
positive_data = torch.utils.data.Subset(dataset, one_indices)


train_loader = DataLoader(positive_data, batch_size=b_size, shuffle=True)

#defining the network
from torch import nn
import torch.nn.functional as F
import torch.optim as optim

class LSTMTagger(nn.Module):

    def __init__(self, embedding_dim, hidden_dim, oligo_set_size):
        super(LSTMTagger, self).__init__()
        self.hidden_dim = hidden_dim
        self.oligo_embeddings = nn.Embedding(oligo_set_size, embedding_dim)

        # The LSTM takes oligonucleotide embeddings as inputs, and outputs hidden states
        # with dimensionality hidden_dim.
        self.lstm = nn.LSTM(embedding_dim, hidden_dim, batch_first=True)

        # The linear layer that maps from hidden state space to tag space
        self.hidden2tag = nn.Linear(hidden_dim, 2)
      
    def forward(self, fragment):
        embeds = self.oligo_embeddings(fragment)
        lstm_out, _ = self.lstm(embeds)
        tag_space = self.hidden2tag(lstm_out)
        tag_scores = F.log_softmax(tag_space, dim=2)
        return tag_scores

#hyper parameters
EMBEDDING_DIM = 50
HIDDEN_DIM = 200
LEARNING_RATE = 0.01
EPOCHS = 300

#train the model
from sklearn.metrics import confusion_matrix

model = LSTMTagger(EMBEDDING_DIM, HIDDEN_DIM, len(dataset.uniq_oligos))
CM=0
loss_function = nn.NLLLoss(reduce=False)
optimizer = optim.SGD(model.parameters(), lr=LEARNING_RATE)

for epoch in range(EPOCHS):
    for j,(fragment_train, targets, weight) in enumerate(train_loader):
        model.zero_grad()
        #calculate output
        tag_scores = model(fragment_train)

        loss = loss_function(tag_scores.view(-1,2), targets.view(-1)) #reshape required due to batch for NLLLoss
        loss = torch.mean(loss * weight.view(-1))
        loss.backward()
        optimizer.step()

        #evaluations 
        preds = torch.argmax(tag_scores, 2) # tag_scores.size() [64, 100, 2]
        CM = confusion_matrix(targets.view(-1), preds.view(-1))
        tn=CM[0][0]
        tp=CM[1][1]
        fp=CM[0][1]
        fn=CM[1][0]
        acc=np.sum(np.diag(CM)/np.sum(CM))
        sensitivity=tp/(tp+fn)
        precision=tp/(tp+fp)
        
#        print('\nAccuracy(mean): %f %%' % (100 * acc))
#        print()
#        print('Confusion Matirx : ')
#        print(CM)
        print('- Sensitivity : ',(tp/(tp+fn))*100, '- Specificity : ',(tn/(tn+fp))*100, '- Precision: ',(tp/(tp+fp))*100)
#        print('- Specificity : ',(tn/(tn+fp))*100)
#        print('- Precision: ',(tp/(tp+fp))*100)
#        print('- NPV: ',(tn/(tn+fn))*100)
#        print('- F1 : ',((2*sensitivity*precision)/(sensitivity+precision))*100)
#        print()
                
#    return acc, CM

#


