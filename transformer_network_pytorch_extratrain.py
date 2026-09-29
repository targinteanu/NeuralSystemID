# %%
import math
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset, ConcatDataset, Dataset
#from sklearn.model_selection import train_test_split
#import matplotlib.pyplot as plt
from myPytorchModels import TimeSeriesConv as myNetMdl
from myPytorchModelTrainer import trainDynsysModel
import pandas as pd
from csv2numpy import prepTimeSeqData

# %%
# Prepare the Data ---------------------------------------------------------------------

seq_len = 128  # sequence length
mdl_Ts = 0.01  # model sample time, s

infofile = "info.csv"
info_df = pd.read_csv(infofile)
feature_names = info_df.iloc[:, 0].values

num_feat = len(feature_names)

# %%
# Initialize the Model, Loss Function, and Optimizer

test_size=0.1
batch_size = 16

groupsize = 15

model = myNetMdl(dim_in=num_feat, dim_out=num_feat, time_len=seq_len, group_size=groupsize, num_groups=5, numGrpUnpaired=2, tuple_size=3)
model.load_state_dict(torch.load("neural_network_pytorch.pth"))

criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

total_params = sum(p.numel() for p in model.parameters() if p.requires_grad)

# %%
# Step 4: Train the Model
optimizer = optim.Adam(model.parameters(), lr=0.001)

# start by performing a few epochs on each set 

train_losses = []
val_losses = []

hzn_len = 16

# %% 
# load final data set 

fs, feature_names, feature_correction, Xs, Ys, X, Y, _, _, _, _, Us, U = prepTimeSeqData(
    seq_len=seq_len, maxNumel=4e9, hzn_len=hzn_len, dt_target=mdl_Ts, 
    filepath="")
Xs = torch.from_numpy(Xs).float()
Ys = torch.from_numpy(Ys).float()
Us = torch.from_numpy(Us).float()
X = torch.from_numpy(X).float()
Y = torch.from_numpy(Y).float()
U = torch.from_numpy(U).float()

# %%

# train / test
# X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size, random_state=42)
train_N = int((1 - test_size) * len(X))
X_train = X[:train_N]
Y_train = Y[:train_N]
U_train = U[:train_N]
X_test = X[train_N:]
Y_test = Y[train_N:]
U_test = U[train_N:]

train_N_s = int((1 - test_size) * len(Xs))
Xs_train = Xs[:train_N_s]
Ys_train = Ys[:train_N_s]
Us_train = Us[:train_N_s]
Xs_test = Xs[train_N_s:]
Ys_test = Ys[train_N_s:]
Us_test = Us[train_N_s:]

# Create TensorDatasets
train_dataset = TensorDataset(X_train, Y_train, U_train)
test_dataset = TensorDataset(X_test, Y_test, U_test)
train_dataset_s = TensorDataset(Xs_train, Ys_train, Us_train)
test_dataset_s = TensorDataset(Xs_test, Ys_test, Us_test)

# Create DataLoaders for batching
train_loader = DataLoader(train_dataset, batch_size, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size, shuffle=False)
#all_loader = DataLoader(TensorDataset(X, Y, U), shuffle=False)
train_loader_s = DataLoader(train_dataset_s, batch_size, shuffle=True)
test_loader_s = DataLoader(test_dataset_s, batch_size, shuffle=False)
#all_loader_s = DataLoader(TensorDataset(Xs, Ys, Us), shuffle=False)

print(f"train dataset size: {len(train_dataset)}, train dataset_s size: {len(train_dataset_s)}")
print(f"test dataset size: {len(test_dataset)}, test dataset_s size: {len(test_dataset_s)}")
print(f"train loader size: {len(train_loader)}, train loader_s size: {len(train_loader_s)}")
print(f"test loader size: {len(test_loader)}, test loader_s size: {len(test_loader_s)}")

print("Total learnable parameters:", total_params)
print("Training data shape:", Y_train.shape, Ys_train.shape)
print("Training data size:", Y_train.numel() + Ys_train.numel())
print("Ratio: ", (Y_train.numel() + Ys_train.numel()) / total_params)

# datasets no longer used once we have loaders, so we can delete to free memory
del train_dataset, test_dataset, train_dataset_s, test_dataset_s
del X_train, Y_train, U_train, X_test, Y_test, U_test, Xs_train, Ys_train, Us_train, Xs_test, Ys_test, Us_test
del X, Y, U, Xs, Ys, Us

# %%
# Step 4A: full training on baseline data 

train_size = len(train_loader.dataset)
steps_per_epoch = math.ceil(train_size / batch_size)
print("Train samples:", train_size)
print("Batch size:", batch_size)
print("Batches/epoch:", steps_per_epoch)

model, train_losses, val_losses = trainDynsysModel(model, optimizer, criterion, train_loader, test_loader, num_epochs=100, allow_early_stopping=True, debugmode=False)

# %%
torch.save(model.state_dict(), "neural_network_pytorch_2.pth")


