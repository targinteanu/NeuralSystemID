# %%
import math
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset, ConcatDataset, Dataset
#from sklearn.model_selection import train_test_split
#import matplotlib.pyplot as plt
from myPytorchModels import TimeSeriesTransformer
from myPytorchModelTrainer import trainDynsysModel
from csv2numpy import prepTimeSeqData

# %%
# Prepare the Data ---------------------------------------------------------------------

seq_len = 64  # sequence length
mdl_Ts = 0.01  # model sample time, s
hzn_len = 16 # samples

fs, feature_names, feature_correction, Xs, Ys, X, Y, _, _, _, _, Us, U = prepTimeSeqData(
    seq_len=seq_len, maxNumel=1e9, hzn_len=hzn_len, dt_target=mdl_Ts, 
    filepath="")
Xs = torch.tensor(Xs, dtype=torch.float32)
Ys = torch.tensor(Ys, dtype=torch.float32)
Us = torch.tensor(Us, dtype=torch.float32)
X = torch.tensor(X, dtype=torch.float32)
Y = torch.tensor(Y, dtype=torch.float32)
U = torch.tensor(U, dtype=torch.float32)

num_feat = len(feature_names)

# %%
# Initialize the Model, Loss Function, and Optimizer

test_size=0.5
batch_size = 32

groupsize = 15

model = TimeSeriesTransformer(dim_in=num_feat, dim_out=num_feat, time_len=seq_len, group_size=groupsize, num_groups=5, numGrpUnpaired=2, tuple_size=3)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

# bias the model frequency prediction to the center of each band 
fcenter = torch.tensor([4,10,27], dtype=torch.float32)
fbias = fcenter.repeat_interleave(groupsize)
fbias = fbias*mdl_Ts*2*math.pi # scale by model sample time and 2pi to convert to radians
with torch.no_grad():
    model.fcoFreq.bias.copy_(fbias)


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
all_loader = DataLoader(TensorDataset(X, Y, U), shuffle=False)
train_loader_s = DataLoader(train_dataset_s, batch_size, shuffle=True)
test_loader_s = DataLoader(test_dataset_s, batch_size, shuffle=False)
all_loader_s = DataLoader(TensorDataset(Xs, Ys, Us), shuffle=False)

print(f"train dataset size: {len(train_dataset)}, train dataset_s size: {len(train_dataset_s)}")
print(f"test dataset size: {len(test_dataset)}, test dataset_s size: {len(test_dataset_s)}")
print(f"train loader size: {len(train_loader)}, train loader_s size: {len(train_loader_s)}")
print(f"test loader size: {len(test_loader)}, test loader_s size: {len(test_loader_s)}")

# %%
hzn_len = 1 # samples
train_loaderlist = [train_loader, train_loader_s]
test_loaderlist = [test_loader, test_loader_s]
# TO DO: think about organizing this as dict instead of list/tuple 

while hzn_len <= 8:

    _, _, _, Xsh, Ysh, Xh, Yh, _, _, _, _, Ush, Uh = prepTimeSeqData(
        seq_len=seq_len, maxNumel=1e9, hzn_len=hzn_len, dt_target=mdl_Ts, 
        filepath="")
    Xsh = torch.tensor(Xsh, dtype=torch.float32)
    Ysh = torch.tensor(Ysh, dtype=torch.float32)
    Ush = torch.tensor(Ush, dtype=torch.float32)
    Xh = torch.tensor(Xh, dtype=torch.float32)
    Yh = torch.tensor(Yh, dtype=torch.float32)
    Uh = torch.tensor(Uh, dtype=torch.float32)

    Xh_train = Xh[:train_N]
    Yh_train = Yh[:train_N]
    Uh_train = Uh[:train_N]
    Xh_test = Xh[train_N:]
    Yh_test = Yh[train_N:]
    Uh_test = Uh[train_N:]

    Xsh_train = Xsh[:train_N_s]
    Ysh_train = Ysh[:train_N_s]
    Ush_train = Ush[:train_N_s]
    Xsh_test = Xsh[train_N_s:]
    Ysh_test = Ysh[train_N_s:]
    Ush_test = Ush[train_N_s:]

    # Create TensorDatasets and loaders
    train_dataset_h = TensorDataset(Xh_train, Yh_train, Uh_train)
    test_dataset_h = TensorDataset(Xh_test, Yh_test, Uh_test)
    train_dataset_sh = TensorDataset(Xsh_train, Ysh_train, Ush_train)
    test_dataset_sh = TensorDataset(Xsh_test, Ysh_test, Ush_test)
    train_loader_h = DataLoader(train_dataset_h, batch_size, shuffle=True)
    test_loader_h = DataLoader(test_dataset_h, batch_size, shuffle=False)
    train_loader_sh = DataLoader(train_dataset_sh, batch_size, shuffle=True)
    test_loader_sh = DataLoader(test_dataset_sh, batch_size, shuffle=False)

    # concat into loader lists 
    train_loaderlist.extend([train_loader_h, train_loader_sh])
    test_loaderlist.extend([test_loader_h, test_loader_sh])

    hzn_len *= 2

train_loaderlist = tuple(train_loaderlist)
test_loaderlist = tuple(test_loaderlist)
print(f"train loader list size: {len(train_loaderlist)}, test loader list size: {len(test_loaderlist)}")

# %%
total_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
print("Total learnable parameters:", total_params)
print("Training data shape:", X_train.shape)
print("Training data size:", X_train.numel())

# %%
# data baseline characteristics as reference for loss 
mean_y = Y_train[:,-1,:].mean(dim=0)
std_y = Y_train[:,-1,:].std(dim=0)
var_y = std_y ** 2
var_per_feat = np.var(Y_train.numpy()[:,-1,:], axis=0)  # redundant?

# %%
# Step 4: Train the Model
optimizer = optim.Adam(model.parameters(), lr=0.001)

# start by performing a few epochs on each set 

#model, trainloss1, valloss1 = trainDynsysModel(model, optimizer, criterion, train_loader, test_loader, num_epochs=5, allow_early_stopping=False)
#model, trainloss2, valloss2 = trainDynsysModel(model, optimizer, criterion, train_loader_s, test_loader_s, num_epochs=5, allow_early_stopping=False)
model, train_losses, val_losses = trainDynsysModel(model, optimizer, criterion, train_loaderlist, test_loaderlist, num_epochs=10, allow_early_stopping=False)
#train_losses = trainloss1 + trainloss2
#val_losses = valloss1 + valloss2

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
# Step 4A: full training on baseline data 

train_size = len(train_loader.dataset)
steps_per_epoch = math.ceil(train_size / batch_size)
print("Train samples:", train_size)
print("Batch size:", batch_size)
print("Batches/epoch:", steps_per_epoch)

model, train_losses, val_losses = trainDynsysModel(model, optimizer, criterion, train_loader, test_loader, num_epochs=100, allow_early_stopping=True)

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
# load or save model
# model.load_state_dict(torch.load("neural_network_pytorch.pth"))
# torch.save(model.state_dict(), "neural_network_pytorch_574259a598db91291cbea59a3d72b242abcdd6b7(4).pth")

# %%
# Step 5A: Evaluate the Model on Test Data

Y_pred = []
Y_test1 = []

model.eval()
with torch.no_grad():
    total_loss = 0
    for x_batch, y_batch, u_batch in test_loader:
        y_pred = model(x_batch, u_batch)
        loss = criterion(y_pred, y_batch)
        total_loss += loss.item() * x_batch.size(0)  # sum up batch loss
        if y_pred.shape[0] == batch_size:
            Y_pred.append(y_pred[:,-1,:])  # only save the final step of the rollout for evaluation
            Y_test1.append(y_batch[:,-1,:])  # only save the final step of the rollout for evaluation

    avg_loss = total_loss / len(test_dataset)
    print(f"Test Loss: {avg_loss:.4f}")

Y_pred_np = np.array(Y_pred)
Y_pred_np = Y_pred_np.reshape(-1, num_feat)
Y_test_np = Y_test[:,-1,:].numpy()
Y_test1_np = np.array(Y_test1)
Y_test1_np = Y_test1_np.reshape(-1, num_feat)

Y_null_all_np = X.numpy()[:, -1, :Y.shape[-1]]
Y_null_test_np = X_test.numpy()[:, -1, :Y.shape[-1]]

MSE_per_feat = np.mean((Y_test1_np - Y_pred_np) ** 2, axis=0)
MSE_per_feat_null = np.mean((Y_test_np - Y_null_test_np) ** 2, axis=0)
feats = np.arange(1, Y.shape[-1]+1)
barwid = .35

"""
plt.figure(figsize=(15,5))
plt.bar(feats - barwid, var_per_feat, width=barwid, label='Output Variance')
plt.bar(feats, MSE_per_feat_null, width=barwid, label='Null MSE')
plt.bar(feats + barwid, MSE_per_feat, width=barwid, label='Test MSE')
plt.xlabel('Output Feature')
plt.xticks(ticks=range(0, len(feats), groupsize), labels=feature_names[::groupsize], rotation=90, ha='right')
plt.ylabel('Value')
plt.title('Output Feature Variance vs Test MSE')
plt.legend()
plt.show()
"""

MSE_per_feat = MSE_per_feat / (np.mean(Y_test1_np**2, axis=0) + np.finfo(float).eps)
MSE_per_feat_null = MSE_per_feat_null / (np.mean(Y_test_np**2, axis=0) + np.finfo(float).eps)

# %%
X_all_np = X.numpy()
Y_all_np = Y[:,-1,:].numpy()
U_all_np = U.numpy()

Y_all_pred = []

model.eval()
with torch.no_grad():
    total_loss = 0
    for x_batch, y_batch, u_batch in all_loader:
        y_pred = model(x_batch, u_batch)
        loss = criterion(y_pred, y_batch)
        total_loss += loss.item() * x_batch.size(0)  # sum up batch loss
        Y_all_pred.append(y_pred[:,-1,:])  # only save the final step of the rollout for evaluation

    avg_loss = total_loss / len(all_loader.dataset)
    print(f"Test+Train Loss: {avg_loss:.4f}")

Y_pred_all_np = np.array(Y_all_pred)
Y_pred_all_np = Y_pred_all_np.reshape(-1, num_feat)
print("Y_pred_all_np shape:", Y_pred_all_np.shape)

# %%
X_all = Xs
Y_all = Ys
U_all = Us

# simulations 
simdur = int(0.2 * fs) # samples 
plotdomain = 1000 * np.array([-1, 1]) + train_N

Ysim = []
i0 = plotdomain[0]
model.eval()
while i0+simdur < plotdomain[1]:
    xi = X_all[i0, :, :].reshape(1,-1,X_all.shape[-1])
    ui = U_all[i0:i0+simdur,0:1,:] # use rollout
    ui = ui.permute(1,0,2)
    with torch.no_grad():
        yi = model(xi, ui)
    Ysim.append(yi[0,:,:].numpy())
    i0 += simdur
    print("Simulating:", (i0-plotdomain[0])/(plotdomain[1]-plotdomain[0]), " complete." )

Ysim = np.concatenate(Ysim, axis=0)
print("Ysim shape:", Ysim.shape)
plotxval = np.arange(len(Ysim)) + plotdomain[0]

# %%
# show several examples 

iMSE = np.argsort(MSE_per_feat)
iVAR = np.argsort(var_per_feat)
iLRN = np.argsort(MSE_per_feat / (var_per_feat + np.finfo(float).eps))
iToPlot = [iMSE[:2], iMSE[-2:], iVAR[:2], iVAR[-2:], iLRN[:2], iLRN[-2:]]
iToPlot = list(set([i for sublist in iToPlot for i in sublist]))
iToPlot = [i-1 for i in iToPlot]  # adjust for zero indexing

"""
plt.figure(figsize=(15,20))
iPlot = 1
for i in iToPlot:
    plt.subplot(len(iToPlot), 1, iPlot)
    plt.plot(Y_all_np[:, i], label='True')
    plt.plot(Y_pred_all_np[:, i], label='Predicted', linestyle='--')
    plt.plot(Y_null_all_np[:, i], label='Null', linestyle=':')
    plt.plot(plotxval, Ysim[:,i], label='Simulated', linestyle='-.')
    plt.xlim(plotdomain)

    # set the y limits to be slightly larger than the min/max of true values in the plotdomain
    y_min = np.min(Y_all_np[plotdomain[0]:plotdomain[1], i])
    y_max = np.max(Y_all_np[plotdomain[0]:plotdomain[1], i])
    y_range = y_max - y_min
    plt.ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)
    
    plt.title(f'Feature {feature_names[i]} - MSE: {MSE_per_feat[i]:.4f}, VAR: {var_per_feat[i]:.4f}')
    plt.legend(loc='upper right')
    iPlot += 1
plt.tight_layout()
plt.show()
"""

# %%
# Step 4B: full training on main data 

train_size = len(train_loader_s.dataset)
steps_per_epoch = math.ceil(train_size / batch_size)
print("Train samples:", train_size)
print("Batch size:", batch_size)
print("Batches/epoch:", steps_per_epoch)

model, train_losses, val_losses = trainDynsysModel(model, optimizer, criterion, train_loader_s, test_loader_s, num_epochs=100, allow_early_stopping=True)

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
# Step 5B: Evaluate the Model on Test Data

Y_pred = []
Y_test1 = []

model.eval()
with torch.no_grad():
    total_loss = 0
    for x_batch, y_batch, u_batch in test_loader_s:
        y_pred = model(x_batch, u_batch)
        loss = criterion(y_pred, y_batch)
        total_loss += loss.item() * x_batch.size(0)  # sum up batch loss
        if y_pred.shape[0] == batch_size:
            Y_pred.append(y_pred[:,-1,:])  # only save the final step of the rollout for evaluation
            Y_test1.append(y_batch[:,-1,:])  # only save the final step of the rollout for evaluation

    avg_loss = total_loss / len(test_dataset_s)
    print(f"Test Loss: {avg_loss:.4f}")

Y_pred_np = np.array(Y_pred)
Y_pred_np = Y_pred_np.reshape(-1, num_feat)
Y_test_np = Ys_test[:,-1,:].numpy()
Y_test1_np = np.array(Y_test1)
Y_test1_np = Y_test1_np.reshape(-1, num_feat)

Y_null_all_np = Xs.numpy()[:, -1, :Ys.shape[-1]]
Y_null_test_np = Xs_test.numpy()[:, -1, :Ys.shape[-1]]

MSE_per_feat = np.mean((Y_test1_np - Y_pred_np) ** 2, axis=0)
MSE_per_feat_null = np.mean((Y_test_np - Y_null_test_np) ** 2, axis=0)
feats = np.arange(1, Ys.shape[-1]+1)
barwid = .35

"""
plt.figure(figsize=(15,5))
plt.bar(feats - barwid, var_per_feat, width=barwid, label='Output Variance')
plt.bar(feats, MSE_per_feat_null, width=barwid, label='Null MSE')
plt.bar(feats + barwid, MSE_per_feat, width=barwid, label='Test MSE')
plt.xlabel('Output Feature')
plt.xticks(ticks=range(0, len(feats), groupsize), labels=feature_names[::groupsize], rotation=90, ha='right')
plt.ylabel('Value')
plt.title('Output Feature Variance vs Test MSE')
plt.legend()
plt.show()
"""

MSE_per_feat = MSE_per_feat / (np.mean(Y_test1_np**2, axis=0) + np.finfo(float).eps)
MSE_per_feat_null = MSE_per_feat_null / (np.mean(Y_test_np**2, axis=0) + np.finfo(float).eps)

# %%
X_all_np = Xs.numpy()
Y_all_np = Ys[:,-1,:].numpy()
U_all_np = Us.numpy()

Y_all_pred = []

model.eval()
with torch.no_grad():
    total_loss = 0
    for x_batch, y_batch, u_batch in all_loader_s:
        y_pred = model(x_batch, u_batch)
        loss = criterion(y_pred, y_batch)
        total_loss += loss.item() * x_batch.size(0)  # sum up batch loss
        Y_all_pred.append(y_pred[:,-1,:])  # only save the final step of the rollout for evaluation

    avg_loss = total_loss / len(all_loader_s.dataset)
    print(f"Test+Train Loss: {avg_loss:.4f}")

Y_pred_all_np = np.array(Y_all_pred)
Y_pred_all_np = Y_pred_all_np.reshape(-1, num_feat)
print("Y_pred_all_np shape:", Y_pred_all_np.shape)

# %%
X_all = Xs
Y_all = Ys
U_all = Us

# simulations 
simdur = int(0.2 * fs) # samples 
plotdomain = 5000 * np.array([-1, 1]) + train_N_s

Ysim = []
i0 = plotdomain[0]
model.eval()
while i0+simdur < plotdomain[1]:
    xi = X_all[i0, :, :].reshape(1,-1,X_all.shape[-1])
    ui = U_all[i0:i0+simdur,0:1,:] # use rollout
    ui = ui.permute(1,0,2)
    with torch.no_grad():
        yi = model(xi, ui)
    Ysim.append(yi[0,:,:].numpy())
    i0 += simdur
    print("Simulating:", (i0-plotdomain[0])/(plotdomain[1]-plotdomain[0]), " complete." )

Ysim = np.concatenate(Ysim, axis=0)
plotxval = np.arange(len(Ysim)) + plotdomain[0]
print("Ysim shape:", Ysim.shape)

# %%
# show several examples 

iMSE = np.argsort(MSE_per_feat)
iVAR = np.argsort(var_per_feat)
iLRN = np.argsort(MSE_per_feat / var_per_feat)
iToPlot = [iMSE[:2], iMSE[-2:], iVAR[:2], iVAR[-2:], iLRN[:2], iLRN[-2:]]
iToPlot = list(set([i for sublist in iToPlot for i in sublist]))

"""
plt.figure(figsize=(15,20))
iPlot = 1
for i in iToPlot:
    plt.subplot(len(iToPlot)+1, 1, iPlot)
    plt.plot(Y_all_np[:, i], label='True')
    plt.plot(Y_pred_all_np[:, i], label='Predicted', linestyle='--')
    plt.plot(Y_null_all_np[:, i], label='Null', linestyle=':')
    plt.plot(plotxval, Ysim[:,i], label='Simulated', linestyle='-.')
    plt.xlim(plotdomain)

    # set the y limits to be slightly larger than the min/max of true values in the plotdomain
    y_min = np.min(Y_all_np[plotdomain[0]:plotdomain[1], i])
    y_max = np.max(Y_all_np[plotdomain[0]:plotdomain[1], i])
    y_range = y_max - y_min
    plt.ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)
    
    plt.title(f'Feature {feature_names[i]} - MSE: {MSE_per_feat[i]:.4f}, VAR: {var_per_feat[i]:.4f}')
    plt.legend(loc='upper right')
    iPlot += 1

plt.subplot(len(iToPlot)+1, 1, iPlot)
plt.plot(U_all_np[:,0,0])
plt.xlim(plotdomain)
plt.title('Stim Input')

plt.tight_layout()
plt.show()
"""

# %%
torch.save(model.state_dict(), "neural_network_pytorch.pth")


