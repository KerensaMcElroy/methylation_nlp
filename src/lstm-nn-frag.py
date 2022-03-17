#importing libraries
import torch
import numpy as np
import matplotlib.pyplot as plt
import sklearn.metrics as metrics

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
        self.meth = self.df[1]

    def get_uniq_oligos(self):
        all_oligos = self.oligos.str.cat(sep=',').split(',')
        oligo_counts = Counter(all_oligos)
        return sorted(oligo_counts, key=oligo_counts.get, reverse=True)
    
    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):

        oligo_list = self.df[0][idx].split(',')
        oligo_indexes = [self.oligo_to_index[o] for o in oligo_list]
        cont_meth = np.array(self.meth[idx].split(',')).astype(np.float) 
        binary_meth = np.where(cont_meth > 0, 1, 0)
        
        return (
            torch.tensor(oligo_indexes),
            torch.tensor(binary_meth)
        )

#create dataloader

meth_file = 'window_mean_HE1_chrm11_iaN7_win10_step10_frag1000.txt'

dataset = DnaMethylation(meth_file)


train_length=int(0.7* len(dataset))
test_length=len(dataset)-train_length
b_size=2

trainset, testset = random_split(dataset, [train_length, test_length])
train_loader = DataLoader(trainset, batch_size=b_size, shuffle=True)
test_loader = DataLoader(testset, batch_size=b_size, shuffle=True)

""" for i, (oligos, meths) in enumerate(train_loader):
     print(oligos.shape, meths.shape)
    break"""

#defining the network
from torch import nn
import torch.nn.functional as F

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
       # print('fragment shape:',fragment.shape)
        embeds = self.oligo_embeddings(fragment)
       # print(len(fragment))
       # print('embed shape:',embeds.shape)
        lstm_out, _ = self.lstm(embeds)
       # print('lstm shape:', lstm_out.shape)
        tag_space = self.hidden2tag(lstm_out)
       # print('tagspace shape:' ,tag_space.shape)
        tag_scores = F.log_softmax(tag_space, dim=2)
       # print('tag_scores shape:',tag_scores.shape)
        #print('tagscores:',tag_scores)
        return tag_scores

#hyper parameters
embedding_dim = 8
hidden_dim = 6
learning_rate = 0.01
epochs = 100

# Model , Optimizer, Loss
model = LSTMTagger(embedding_dim, hidden_dim, len(dataset.uniq_oligos))
optimizer = torch.optim.SGD(model.parameters(),lr=learning_rate)
loss_function = nn.NLLLoss()

#train the model

losses = []
accur = []
truepr = []

for i in range(epochs):
  for j,(fragment_train,meth_train) in enumerate(train_loader):
    model.zero_grad()
    #calculate output
    meth_predict = model(fragment_train)
#    print('meth_predict shape:', meth_predict.shape)
    loss = loss_function(meth_predict.view(-1,2), meth_train.view(-1))

    
    #backprop
    loss.backward()
    optimizer.step()

    #evaluation
    #acc = (meth_predict.reshape(-1).detach().numpy().round() == meth_train.mean()
    positives = np.sum(meth_train.reshape(-1).detach().numpy() == 1)
    value, prediction = torch.max(meth_train, 2)
    print('pred shape', prediction.shape)
    print(prediction)
    predict_neg = meth_predict[:,:,1].reshape(-1).detach().numpy().round() == 0. 
   # print(meth_predict[:,:,1].reshape(-1).detach().numpy()) 
   # tp = np.sum(( & (meth_train.reshape(-1).detach().numpy() == 1))
  #  fn = np.sum((meth_predict[:,:,1].reshape(-1).detach().numpy().round() == 0) & (meth_train.reshape(-1).detach().numpy() == 1))
   # tpr = tp/(tp+fn)

    if i%50 == 0:
        losses.append(loss.detach().numpy())
       # accur.append(acc)
       # truepr.append(tpr)
       # print("epoch: {}\tloss: {}\ttpr: {}".format(i, loss, tpr))# {}\t accuracy : {}\t tpr: {}".format(i,loss,acc, tpr))

#plotting results
#plt.plot(losses)
#plt.title('Loss vs Epochs')
#plt.xlabel('Epochs')
#plt.ylabel('loss')

#plt.savefig("losses.png")

#plt.plot(accur)
#plt.title('Accuracy vs Epochs')
#plt.xlabel('Epochs')
#plt.ylabel('accuracy')

#plt.savefig("accuracy.png")

#plt.plot(truepr)
#plt.title('True Positive Rate vs Epochs')
#plt.xlabel('Epochs')
#plt.ylabel('TPR')

#plt.savefig("tpr.png")

