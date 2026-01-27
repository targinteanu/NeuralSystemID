import math
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import pandas as pd
import matplotlib.pyplot as plt
import bisect
from myPytorchModels import TimeSeriesTransformer

# set params -------------------------------------------------------------------------------------
hzn = .05 # EVALUATION sample time, s
groupsize=16
numgroups=4
f = np.array([4,10,27,70]) # freq band center freqs
netfile = "neural_network_pytorch_1639e948afec6ecef1e37515d791e25403ee7dd4.pth"
dt_target = 0.005 # model sample time, s
seq_len = 64 # model transformer samples
hzn_len = math.ceil(hzn / dt_target)  # horizon as multiple of MODEL Ts, NOT data Ts 

# Prepare the Data ---------------------------------------------------------------------

# ------------------------
# Load info
# ------------------------
info_df = pd.read_csv("info.csv")
fs = info_df.iloc[0, 5]                  # sampling frequency (Hz)
print("Sampling frequency (Hz):", fs)
feature_names = info_df.iloc[:, 0].values
#feature_correction = info_df.iloc[:, 3].values

# ------------------------
# Load raw data
# ------------------------
baseline_df = pd.read_csv("baselinedataraw.csv")
baseline_time = baseline_df.iloc[:, 0].values            # time column (seconds)
baseline_data = baseline_df.iloc[:, 1:].values           # data columns

iBLstart = math.floor(0.995*baseline_data.shape[0])
baseline_time = baseline_time[iBLstart:]
baseline_data = baseline_data[iBLstart:, :]

df = pd.read_csv("dataraw.csv")
time = df.iloc[:, 0].values            # time column (seconds)
data = df.iloc[:, 1:].values           # data columns

# ------------------------
# Determine outliers
# ------------------------
threshsd = 3 # standard deviations 
threshprop = .5 # proportion of features
BLmean = np.mean(baseline_data, axis=0)
BLstd = np.std(baseline_data, axis=0)
BLisout = np.abs(baseline_data - BLmean) > (threshsd * BLstd)
BLisnoise = np.sum(BLisout, axis=1) > (threshprop * baseline_data.shape[1])
isout = np.abs(data - BLmean) > (threshsd * BLstd)
isnoise = np.sum(isout, axis=1) > (threshprop * data.shape[1])

# ------------------------
# Load processed data
# ------------------------
baseline_df = pd.read_csv("baselinedata.csv")
baseline_time = baseline_df.iloc[:, 0].values            # time column (seconds)
baseline_data = baseline_df.iloc[:, 1:].values           # data columns

baseline_time = baseline_time[iBLstart:]
baseline_data = baseline_data[iBLstart:, :]

df = pd.read_csv("data.csv")
time = df.iloc[:, 0].values            # time column (seconds)
data = df.iloc[:, 1:].values           # data columns

# ------------------------
# Load events
# ------------------------
events_df = pd.read_csv("events.csv")
event_times = events_df.iloc[:, 0].values  # assume first column is event time in seconds
event_times = np.sort(event_times)         # ensure sorted

# For fast lookup using binary search
def count_events_in_window(t, window=0.2):
    """
    Count how many event_times fall in (t - window, t].
    Uses bisect for O(log n) search.
    """
    left = bisect.bisect_right(event_times, t - window)
    right = bisect.bisect_right(event_times, t)
    return right - left

# ------------------------
# Compute event count for each row
# ------------------------
event_counts = np.array([count_events_in_window(t, 1.2/fs) for t in time])
event_counts = event_counts.reshape(-1, 1)

# Append event_counts as an additional input feature
data_aug = np.hstack([data, event_counts])
# Now each input row has: [original data..., event_count]

# ------------------------
# Build input-output pairs 
# ------------------------
dt_tol = 0.15 * dt_target
drow_target = int(dt_target * fs)  # number of rows 

X_list = []
Y_list = []

# create sliding windows
def create_windows(data, seq_len=128, horizon=1):
    X, Y = [], []
    for i in range(len(data) - seq_len - horizon + 1):
        X.append(data[i:i+seq_len])
        Y.append(data[i+seq_len + horizon - 1])
    return np.array(X), np.array(Y)

# --- baseline data ---
#Nbl = len(baseline_df)
Nbl = len(baseline_data)
for iStart in range(drow_target):
    inputs = []
    for i in range(iStart, (Nbl - drow_target), drow_target):
        dt = baseline_time[i+drow_target] - baseline_time[i]
        if (abs(dt - dt_target) <= dt_tol) and (not BLisnoise[i]):
            inputs.append(baseline_data[i, :]) 
        else:
            if len(inputs) > seq_len+hzn_len:
                x, y = create_windows(inputs, seq_len, hzn_len)
                if x is not None:
                    X_list.append(x)
                    Y_list.append(y)
            inputs = []
    # catch trailing segment
    if len(inputs) > seq_len+hzn_len:
        x, y = create_windows(inputs, seq_len, hzn_len)
        if x is not None:
            X_list.append(x)
            Y_list.append(y)

X = np.concatenate(X_list, axis=0)
Y = np.concatenate(Y_list, axis=0)

X = torch.tensor(X, dtype=torch.float32)
Y = torch.tensor(Y, dtype=torch.float32)

print("Pairs created:", len(X))
print("Input shape :", X.shape)   
print("Output shape:", Y.shape)

num_feat = X.shape[-1]

# -----------------------------------------------------------------------------------------------

# define mdl struct ====================================================================
model = TimeSeriesTransformer(dim_in=num_feat, dim_out=num_feat, time_len=seq_len, group_size=groupsize, num_groups=numgroups)
model.load_state_dict(torch.load(netfile))

# simulations -----------------------------------------------------------------------------------
Ysim = []
model.eval()
progtick = 0.05
progcur = 0.0
prognext = progtick
print("Simulating...")
for i0 in range(len(X) - hzn_len):
    xi = X[i0, :, :].reshape(1,-1,num_feat)
    for i in range(hzn_len):
        with torch.no_grad():
            yi = model(xi).numpy().flatten()
        # prepare next input
        if i < hzn_len - 1:
            #event_count_next = X_all_np[i0 + i + 1, -1]  # keep using original event count
            #xi = torch.tensor(np.hstack([yi, event_count_next]).reshape(1, -1), dtype=torch.float32)
            xi = torch.tensor(np.vstack([xi[0, 1:, :], yi]).reshape(1,-1,X.shape[-1]), dtype=torch.float32)
    Ysim.append(yi)
    progcur += 1.0/(len(X) - hzn_len)
    if progcur >= prognext:
        print(f"  {int(progcur*100)}%")
        prognext += progtick

# evaluate -------------------------------------------------------------------------------------

# overall 
Ysim = np.array(Ysim)
Ytrue = Y[hzn_len:, :].numpy()
print("Ysim shape :", Ysim.shape)
print("Ytrue shape:", Ytrue.shape)
mse = np.mean((Ysim - Ytrue)**2, axis=0) / np.mean((Ytrue)**2, axis=0)
rho = np.array([ np.corrcoef(Ysim[:,i], Ytrue[:,i])[0,1] for i in range(Ytrue.shape[1]) ])
print("Mean MSE:", np.mean(mse))
print("Mean correlation:", np.mean(rho))
# show a bar plot of mse per feature
barwid = .35
barx = np.arange(len(mse))
plt.figure()
plt.bar(barx-barwid, 1-mse, width=barwid, label="1-MSE")
plt.bar(barx, rho, width=barwid, label="correlation")
plt.legend()
plt.xlabel("Feature index")
plt.xticks(ticks=range(0, len(mse), groupsize), labels=feature_names[::groupsize], rotation=90, ha='right')
plt.ylabel("accuracy")
plt.title("Accuracy per Feature")
plt.grid(axis='y')
plt.show()

# a few example channels 
# get index of channels sorted by mse
sorted_indices = np.argsort(rho-mse)  # ascending order
# get first, middle, and end of sorted_indices
example_indices = [sorted_indices[0], sorted_indices[len(sorted_indices)//4], sorted_indices[len(sorted_indices)//2], sorted_indices[3*len(sorted_indices)//4], sorted_indices[-1]]
# Create a single figure with vertically stacked subplots
fig, axes = plt.subplots(len(example_indices), 1, sharex=True, figsize=(10, 8))
for ax, idx in zip(axes, example_indices):
    ax.plot(Ytrue[:, idx], label="True")
    ax.plot(Ysim[:, idx], label="Simulated", alpha=0.7)
    ax.set_ylabel("Feature Value")
    ax.set_title(f"Feature: {feature_names[idx]} (MSE: {mse[idx]:.4f}; corr: {rho[idx]:.4f})")
    ax.legend()
    ax.grid(axis='both')
axes[-1].set_xlabel("Sample")  # Set x-label only on the last subplot
plt.tight_layout()
plt.show()

# reconstruct raw signal from filtered bands 
Ysim_grouped = Ysim[:, :numgroups * groupsize] # * feature_correction[:numgroups * groupsize]
Ytrue_grouped = Ytrue[:, :numgroups * groupsize] # * feature_correction[:numgroups * groupsize]
mse_grouped = mse[:numgroups * groupsize].reshape(groupsize, numgroups)
"""
Yimagsim_grouped = Ysim[:, (numgroups*groupsize):(2*numgroups*groupsize)].reshape(-1, groupsize, numgroups)
Yimagtrue_grouped = Ytrue[:, (numgroups*groupsize):(2*numgroups*groupsize)].reshape(-1, groupsize, numgroups)
print("Yimagsim_grouped shape : ", Yimagsim_grouped.shape)
print("Yimagtrue_grouped shape: ", Yimagtrue_grouped.shape)
Psim = math.log10( Ysim_grouped**2 + Yimagsim_grouped**2 )
Ptrue = math.log10( Ytrue_grouped**2 + Yimagtrue_grouped**2 )
"""
pinksim = Ysim[:, 2*(numgroups * groupsize):]
pinktrue = Ytrue[:, 2*(numgroups * groupsize):]
f = np.tile(f, (groupsize,1)).T.flatten()
f = np.column_stack((np.log10(f), np.ones(len(f)))).T
Ppinksim = np.power(10, pinksim @ f) ** .5
Ppinktrue = np.power(10, pinktrue @ f) ** .5
Ysim_grouped = Ysim_grouped * Ppinksim
Ytrue_grouped = Ytrue_grouped * Ppinktrue
Ysim_grouped = Ysim_grouped.reshape(-1, numgroups, groupsize)
Ytrue_grouped = Ytrue_grouped.reshape(-1, numgroups, groupsize)
Ysim_recon = np.sum(Ysim_grouped, axis=-2)
Ytrue_recon = np.sum(Ytrue_grouped, axis=-2)
mse_recon = np.mean((Ysim_recon - Ytrue_recon)**2, axis=0) / np.mean((Ytrue_recon)**2, axis=0)
rho_recon = np.array([ np.corrcoef(Ysim_recon[:,i], Ytrue_recon[:,i])[0,1] for i in range(Ytrue_recon.shape[1]) ])

# show a bar plot of mse per feature
barx = np.arange(len(mse_recon))
plt.figure()
plt.bar(barx-barwid, 1-mse_recon, width=barwid, label="1-MSE")
plt.bar(barx, rho_recon, width=barwid, label="correlation")
plt.legend()
plt.xlabel("Channel index")
plt.ylabel("accuracy")
plt.title("Accuracy per Channel")
plt.grid(axis='y')
plt.show()

# a few example channels 
# get index of channels sorted by mse
sorted_indices = np.argsort(rho_recon-mse_recon)  # ascending order
# get first, middle, and end of sorted_indices
example_indices = [sorted_indices[0], sorted_indices[len(sorted_indices)//4], sorted_indices[len(sorted_indices)//2], sorted_indices[3*len(sorted_indices)//4], sorted_indices[-1]]
# Create a single figure with vertically stacked subplots
fig, axes = plt.subplots(len(example_indices), 1, sharex=True, figsize=(10, 8))
for ax, idx in zip(axes, example_indices):
    ax.plot(Ytrue_recon[:, idx], label="True")
    ax.plot(Ysim_recon[:, idx], label="Simulated", alpha=0.7)
    ax.set_ylabel("Channel Value")
    ax.set_title(f"Channel: {idx} (MSE: {mse_recon[idx]:.4f}; corr: {rho_recon[idx]:.4f})")
    ax.legend()
    ax.grid(axis='both')
axes[-1].set_xlabel("Sample")  # Set x-label only on the last subplot
plt.tight_layout()
plt.show()