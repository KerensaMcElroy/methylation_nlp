import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np
from torch.utils.data import DataLoader, random_split, Dataset
import time
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

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
        cont_meth = np.loadtxt(src_file, max_rows=num_rows,
            usecols=2, delimiter=" ", skiprows=0,
            dtype=np.float32)
        self.m_data = np.where(cont_meth > 0, 1, 0)
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
        meth = [self.m_data[idx]]
        start = self.s_data[idx]
        end = self.e_data[idx]
    
        return (oligo, meth, start, end)


def prepare_sequence(seq, to_ix):
    idxs = [to_ix[s] for s in seq]
    idxs = torch.tensor(idxs, dtype=torch.long)
    return idxs

dataset = MethDataset('/datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/users/mce03b/methylation_nlp/analysis/window_mean_HE1_chrm2_iaN7_win50_step25_noN_10000.txt')
trainset, valset = random_split(dataset, [round(len(dataset)*0.9), round(len(dataset)*0.1)])
train_loader = DataLoader(trainset, batch_size=1, shuffle=True, num_workers=1)
val_loader = DataLoader(valset, batch_size=1, shuffle=True, num_workers=1)



nuc_to_ix = {"A":0,"T":1,"G":2,"C":3}

class LSTMTagger(nn.Module):

    def __init__(self, nuc_embedding_dim, 
                 nuc_hidden_dim, nn_dim, nuc_vocab_size):

        super(LSTMTagger, self).__init__()
        self.nuc_hidden_dim = nuc_hidden_dim
        self.nuc_embeddings = nn.Embedding(nuc_vocab_size, nuc_embedding_dim)
        self.nuc_hidden = self.init_hidden(self.nuc_hidden_dim)

        # The LSTM takes nucleotide embeddings as inputs, and outputs hidden states
        # with dimensionality nuc_hidden_dim.

        self.nuc_lstm = nn.LSTM(nuc_embedding_dim, nuc_hidden_dim)

        # create fully connected neural network to get from nuc to oligo
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(nuc_hidden_dim, nn_dim),
            nn.ReLU(),
            nn.Linear(nn_dim, nn_dim),
            nn.ReLU(),
            nn.Linear(nn_dim,2),
        )

    def init_hidden(self, size):
        # The axes semantics are (num_layers, minibatch_size, hidden_dim)
        return (torch.zeros(1, 1, size),
                torch.zeros(1, 1, size))

    def forward(self, nuc_sequence):

        nuc_embeds = self.nuc_embeddings(nuc_sequence)
        nuc_lstm_out, self.nuc_hidden = self.nuc_lstm(nuc_embeds.view(len(nuc_sequence), 1, -1), self.nuc_hidden)

        x = self.flatten(nuc_lstm_out[-1])
        tag_space = self.linear_relu_stack(x)
        tag_scores = F.log_softmax(tag_space, dim=1)
        return tag_scores

NN_DIM = 4
NUC_EMBEDDING_DIM = 4
NUC_HIDDEN_DIM = 50
model = LSTMTagger(NUC_EMBEDDING_DIM, NUC_HIDDEN_DIM, NN_DIM,
                   len(nuc_to_ix))

weights = torch.tensor([.001,.999])
loss_function = nn.CrossEntropyLoss(weight=weights)
optimizer = optim.SGD(model.parameters(), lr=0.005)

#with torch.no_grad():
#    for i, (oligo, meth, start, end) in enumerate(train_loader):
#        nuc_sequence = prepare_sequence(oligo, nuc_to_ix)
#        print(nuc_sequence)
#        score = model(nuc_sequence)
#        print(score)

# Keep track of losses for plotting
current_loss = 0
all_losses = []
print_every = 1000
plot_every = 100
num_train = len(trainset)

def timeSince(since):
    now = time.time()
    s = now - since
    m = math.floor(s / 60)
    s -= m * 60
    return '%dm %ds' % (m, s)

start = time.time()

#should we instead be training with random length sequence divided into much smaller oligos?
def train(target_meth, oligo_tensor):
    model.zero_grad()
    model.nuc_hidden = model.init_hidden(model.nuc_hidden_dim)
    predict_meth= model(oligo_tensor)
    loss = loss_function(predict_meth, target_meth)
    loss.backward(retain_graph=True)
    optimizer.step()
    
    return predict_meth, loss.item()

    output, loss = train(category_tensor, line_tensor)
    current_loss += loss
for i, (oligo, meth, start, end) in enumerate(train_loader):
    oligo_tensor = prepare_sequence(oligo[0], nuc_to_ix)
    target_meth = meth[0]
    target_meth = torch.tensor(target_meth)
    predict_meth, loss = train(target_meth, oligo_tensor)
    current_loss += loss

 #   if i % plot_every == 0:
 #       mean_loss = current_loss/plot_every
 #       print(mean_loss)
 #       current_loss = 0
    if target_meth.item()==1:
         print('%d %d%% (%s) %.4f %s / %s' % (i, i / num_train * 100, timeSince(start), loss, predict_meth, target_meth))
 #  if i % plot_every == 0:
  #      mean_loss = current_loss/plot_every
  #      print(mean_loss)
    # Add current loss avg to list of losses
  #  if i % plot_every == 0:
  ##      all_losses.append(current_loss / plot_every)
   #     current_loss = 0

plt.figure()
plt.plot(all_losses)

#with torch.no_grad():
#  for i, (oligo, meth, start, end) in enumerate(val_loader):#        
#    nuc_sequence = prepare_sequence(oligo[0], nuc_to_ix)
#    score = model(nuc_sequence)
#    print(score, meth[0])