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
# simpler model as baseline for comparison 
def trainAR(X, Y):
    Mar = []
    for f in range(Y.shape[-1]):
        x = X[:,:,f]
        y = Y[:,0,f]
        A = np.linalg.lstsq(x, y, rcond=None)[0]
        Mar.append(A)
    Mar = np.stack(Mar, axis=1)
    return Mar

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
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

# bias the model frequency prediction to the center of each band 
fcenter = torch.tensor([4,10,27], dtype=torch.float32)
fbias = fcenter.repeat_interleave(groupsize)
fbias = fbias*mdl_Ts*2*math.pi # scale by model sample time and 2pi to convert to radians
fbias = torch.atanh(fbias / (math.pi + 1e-5)) # apply inverse tanh to get bias in pre-tanh space
with torch.no_grad():
    model.fcoFreq.bias.copy_(fbias)

total_params = sum(p.numel() for p in model.parameters() if p.requires_grad)

# %%
# Step 4: Train the Model
optimizer = optim.Adam(model.parameters(), lr=0.001)

# start by performing a few epochs on each set 

train_losses = []
val_losses = []

# del the _s loaders to free memory
# del train_dataset_s, test_dataset_s, train_loader_s, test_loader_s, Xs_train, Ys_train, Us_train, Xs_test, Ys_test, Us_test

# now the loop
hzn_len = 1
while hzn_len < 16:
    _, _, _, Xsh, Ysh, Xh, Yh, _, _, _, _, Ush, Uh = prepTimeSeqData(
        seq_len=seq_len, maxNumel=4e9, hzn_len=hzn_len, dt_target=mdl_Ts, 
        filepath="")
    Xsh = torch.from_numpy(Xsh).float()
    Ysh = torch.from_numpy(Ysh).float()
    Ush = torch.from_numpy(Ush).float()
    Xh = torch.from_numpy(Xh).float()
    Yh = torch.from_numpy(Yh).float()
    Uh = torch.from_numpy(Uh).float()

    train_N = int((1 - test_size) * len(Xh))
    Xh_train = Xh[:train_N]
    Yh_train = Yh[:train_N]
    Uh_train = Uh[:train_N]
    Xh_test = Xh[train_N:]
    Yh_test = Yh[train_N:]
    Uh_test = Uh[train_N:]

    train_N_s = int((1 - test_size) * len(Xsh))
    Xsh_train = Xsh[:train_N_s]
    Ysh_train = Ysh[:train_N_s]
    Ush_train = Ush[:train_N_s]
    Xsh_test = Xsh[train_N_s:]
    Ysh_test = Ysh[train_N_s:]
    Ush_test = Ush[train_N_s:]

    # train simpler baseline model for comparison 
    if hzn_len == 1:
        Mar = trainAR(Xh_train, Yh_train)
        np.save("neural_network_AR.npy", Mar)
        del Mar

    # Create TensorDatasets and loaders
    train_dataset_h = TensorDataset(Xh_train, Yh_train, Uh_train)
    test_dataset_h = TensorDataset(Xh_test, Yh_test, Uh_test)
    train_dataset_sh = TensorDataset(Xsh_train, Ysh_train, Ush_train)
    test_dataset_sh = TensorDataset(Xsh_test, Ysh_test, Ush_test)
    train_loader_h = DataLoader(train_dataset_h, batch_size, shuffle=True)
    test_loader_h = DataLoader(test_dataset_h, batch_size, shuffle=False)
    train_loader_sh = DataLoader(train_dataset_sh, batch_size, shuffle=True)
    test_loader_sh = DataLoader(test_dataset_sh, batch_size, shuffle=False)

    print("Total learnable parameters:", total_params)
    print("Training data shape:", Yh_train.shape, Ysh_train.shape)
    print("Training data size:", Yh_train.numel() + Ysh_train.numel())
    print("Ratio: ", (Yh_train.numel() + Ysh_train.numel()) / total_params)

    # train
    model, tl, vl = trainDynsysModel(model, optimizer, criterion, (train_loader_h, train_loader_sh), (test_loader_h, test_loader_sh), num_epochs=5, allow_early_stopping=False)
    train_losses.extend(tl)
    val_losses.extend(vl)

    # del to free memory
    del Xsh, Ysh, Ush, Xh, Yh, Uh, Xh_train, Yh_train, Uh_train, Xh_test, Yh_test, Uh_test, Xsh_train, Ysh_train, Ush_train, Xsh_test, Ysh_test, Ush_test
    del train_dataset_h, test_dataset_h, train_dataset_sh, test_dataset_sh, train_loader_h, test_loader_h, train_loader_sh, test_loader_sh

    hzn_len *= 2

"""
# After loop: plot train/val loss to inspect convergence
plt.plot(train_losses, label='train_loss')
plt.plot(val_losses, label='val_loss')
plt.xlabel('Epoch')
plt.ylabel('MSE')
plt.legend()
plt.show()
"""

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

"""
# After loop: plot train/val loss to inspect convergence
plt.plot(train_losses, label='train_loss')
plt.plot(val_losses, label='val_loss')
plt.xlabel('Epoch')
plt.ylabel('MSE')
plt.legend()
plt.show()
"""

# %%
# Step 4B: full training on main data 
"""
train_size = len(train_loader_s.dataset)
steps_per_epoch = math.ceil(train_size / batch_size)
print("Train samples:", train_size)
print("Batch size:", batch_size)
print("Batches/epoch:", steps_per_epoch)

model, train_losses, val_losses = trainDynsysModel(model, optimizer, criterion, train_loader_s, test_loader_s, num_epochs=100, allow_early_stopping=True, debugmode=False)
"""
"""
# After loop: plot train/val loss to inspect convergence
plt.plot(train_losses, label='train_loss')
plt.plot(val_losses, label='val_loss')
plt.xlabel('Epoch')
plt.ylabel('MSE')
plt.legend()
plt.show()
"""

# %%
torch.save(model.state_dict(), "neural_network_pytorch.pth")


