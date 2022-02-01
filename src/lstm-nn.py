import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

torch.manual_seed(1)
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f'Using {device} device')

def prepare_sequence(seq, to_ix):
    idxs = [to_ix[s] for s in seq]
    idxs = torch.tensor(idxs, dtype=torch.long)
    return idxs

training_data = [
    # Tags are now continuous on [0,1]
    (["AATCGAT","GGCTGTG","ATGCTGA","GGGGCGG"], [.2,.94,.05,.75]),
    (["GGCGGCC","AGTGGCG","GGGTATA","ATGCCGT"], [.82,.64,.98,.32])
]


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
            nn.Linear(nn_dim,1),
        )

    def init_hidden(self, size):
        # The axes semantics are (num_layers, minibatch_size, hidden_dim)
        return (torch.zeros(1, 1, size),
                torch.zeros(1, 1, size))

    def forward(self, nuc_sequence):

        nuc_embeds = self.nuc_embeddings(nuc_sequence)
        nuc_lstm_out, self.nuc_hidden = self.nuc_lstm(nuc_embeds.view(len(nuc_sequence), 1, -1), self.nuc_hidden)

        x = self.flatten(nuc_lstm_out[-1])
        score = self.linear_relu_stack(x)
        return score

NN_DIM = 6
CHAR_EMBEDDING_DIM = 6
CHAR_HIDDEN_DIM = 6

model = LSTMTagger(CHAR_EMBEDDING_DIM, CHAR_HIDDEN_DIM, NN_DIM,
                   len(nuc_to_ix))

loss_function = nn.MSELoss()
optimizer = optim.SGD(model.parameters(), lr=0.1)

with torch.no_grad():
    for oligo in training_data[0][0]:
        nuc_sequence = prepare_sequence(oligo, nuc_to_ix)
        print(nuc_sequence)
        score = model(nuc_sequence)
        print(score)

for epoch in range(300):
    for sentence, targets in training_data:
        model.zero_grad()

        for index, oligo in enumerate(sentence):
            # nuc hidden init by oligo
            model.nuc_hidden = model.init_hidden(model.nuc_hidden_dim)

            nuc_sequence = prepare_sequence(oligo, nuc_to_ix)
            target = targets[index]
            target = torch.tensor(target)
            score = model(nuc_sequence)

            loss = loss_function(score, target)
            loss.backward(retain_graph=True)

        optimizer.step()


with torch.no_grad():
    for oligo in training_data[0][0]:
        nuc_sequence = prepare_sequence(oligo, nuc_to_ix)
        score = model(nuc_sequence)
        print(score)