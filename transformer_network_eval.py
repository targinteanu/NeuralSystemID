import math
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import pandas as pd
import matplotlib.pyplot as plt
import bisect

# set params -------------------------------------------------------------------------------------
hzn = 0.05 # EVALUATION sample time, s
groupsize=7
numgroups=6
netfile = "neural_network_pytorch_5bf2fceb0050b6a9d3ca5c5321a891eb0fa26a9c.pth"
dt_target = 0.005 # model sample time, s
seq_len = 32 # model transformer samples
hzn_len = math.ceil(hzn / dt_target)  # horizon as multiple of MODEL Ts, NOT data Ts 

# Prepare the Data ---------------------------------------------------------------------

# ------------------------
# Load info
# ------------------------
info_df = pd.read_csv("info.csv")
fs = info_df.iloc[0, 5]                  # sampling frequency (Hz)
feature_names = info_df.iloc[:, 0].values

# ------------------------
# Load raw data
# ------------------------
baseline_df = pd.read_csv("baselinedataraw.csv")
baseline_time = baseline_df.iloc[:, 0].values            # time column (seconds)
baseline_data = baseline_df.iloc[:, 1:].values           # data columns

iBLstart = math.ceil(0.9*baseline_data.shape[0])
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

class SinusoidalPositionalEncoding(nn.Module):
    def __init__(self, dim_model, max_len=5000):
        super().__init__()

        # Create matrix of shape (max_len, dim_model)
        pe = torch.zeros(max_len, dim_model)
        position = torch.arange(0, max_len).unsqueeze(1)  # (max_len, 1)

        # Divide by log-based frequencies
        div_term = torch.exp(torch.arange(0, dim_model, 2) * (-math.log(10000.0) / dim_model))

        # Apply sin to even indices, cos to odd
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)

        # Register as buffer so it's saved with model but not trainable
        self.register_buffer('pe', pe.unsqueeze(0))  # (1, max_len, dim_model)

    def forward(self, x):
        """
        x: (B, T, dim_model)
        returns: x + positional_encoding[:, :T, :]
        """
        T = x.size(1)
        return x + self.pe[:, :T, :]

class TimeSeriesTransformer(nn.Module):
    def __init__(self, dim_in, dim_out, pos_len=512, group_size=7, num_groups=7):
        super().__init__()

        # transformer properties
        dim_model=32
        num_heads=4
        num_layers=8
        dim_ff=128

        # time-independent preprocessing features -------------------------------------------
        
        C1 = 32
        C2 = 32

        self.num_groups = num_groups
        self.group_size = group_size
        self.num_pairs = num_groups * group_size  # 49
        self.used_for_pairing = self.num_pairs * 2
        self.leftover_dim = dim_in - self.used_for_pairing
        self.pair_output = 8

        # Stage 1: pairwise linear 
        self.pair_fc1 = nn.Linear(2, C1)  # like conv kernel_size=2
        self.pair_fc2 = nn.Linear(C1, self.pair_output)  # like conv kernel_size=1

        # Stage 2: linear over flattened group
        self.group_fcA = nn.Linear(group_size * self.pair_output, C2)
        self.group_fcB = nn.Linear(num_groups * self.pair_output, C2)
        # MLP input = 56 + 56 + leftover(3) = 115
        mlp_in = (num_groups * C2) + (group_size * C2) + self.leftover_dim

        # stage 3: a few MLP layers to get to dim_model
        self.fc1 = nn.Linear(mlp_in, 256)
        self.fc2 = nn.Linear(256, 256)
        self.fc3 = nn.Linear(256, 128)
        self.fc4 = nn.Linear(128, 64)
        self.fc5 = nn.Linear(64, dim_model)

        # ---------------------------------------------------------------------------------------

        # Positional embedding
        #self.pos_emb = nn.Parameter(torch.randn(1, pos_len, dim_model))
        self.pos_emb = SinusoidalPositionalEncoding(dim_model, max_len=pos_len)

        # Transformer Encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=dim_model,
            nhead=num_heads,
            dim_feedforward=dim_ff,
            activation='gelu',
            batch_first=True  # lets inputs be (B, T, D)
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)

        # Output head for next-step prediction
        self.fco1 = nn.Linear(dim_model, 64)
        self.fco2 = nn.Linear(64, 64)
        self.fco3 = nn.Linear(64, 64)
        self.fco4 = nn.Linear(64, dim_out)

    def forward(self, x):
        """
        x: (B, T, dim_in)
        """
        if x.ndim != 3:
            raise ValueError(f"Expected input ndim=3, got {x.ndim}")
        T = x.size(1)
        #if T > self.pos_emb.size(1):
        #    raise RuntimeError(f"Sequence length T={T} exceeds pos_len={self.pos_emb.size(1)}. "
        #                       "Either increase pos_len or ensure inputs have smaller T.")
        #x = self.input_proj(x) + self.pos_emb[:, :T, :]
        B = x.size(0)

        x_left = x[:, :, self.used_for_pairing:]  # (B,T,3)
        x_used = x[:, :, :self.used_for_pairing]  # (B, 100)
        x_pairs = torch.stack( 
            (
                x_used[:, :, :self.num_pairs],            # x[0..49]
                x_used[:, :, self.num_pairs:],             # x[49..100]
            ),
            dim=3
        )  # (B,T,49,2)

        # Stage 1
        p = F.gelu(self.pair_fc1(x_pairs))     # (B,T,49,8)
        p = F.gelu(self.pair_fc2(p))           # (B,T,49,2)

        # Stage 2A: groups
        p_groups = p.view(B, T, self.num_groups, self.group_size * self.pair_output)  # (B,T,7,14)
        a = F.gelu(self.group_fcA(p_groups))                      # (B,T,7,8)
        a_flat = a.view(B, T, -1)                                 # (B,T,56)

        # Stage 2B: threads
        p_threads = p.view(B, T, self.num_groups, self.group_size, self.pair_output)   # (B,T,7,7,2)
        p_threads = p_threads.permute(0,1,3,2,4).contiguous()           # (B,T,7,7,2)
        p_threads = p_threads.view(B, T, self.group_size, -1)            # (B,T,7,14)
        b = F.gelu(self.group_fcB(p_threads))                         # (B,T,7,8)
        b_flat = b.view(B, T, -1)                                        # (B,T,56)

        # Stage 3 MLP
        h = torch.cat([a_flat, b_flat, x_left], dim=2)        # (B,T,115)
        h = F.gelu(self.fc1(h))
        h = F.gelu(self.fc2(h))
        h = F.gelu(self.fc3(h))
        h = F.gelu(self.fc4(h))
        h = F.gelu(self.fc5(h))

        # Transformer --------------------------------------------------------------------

        h = self.pos_emb(h)
        z = self.encoder(h)

        # Output head
        y = z[:, -1]  # (B, dim_model)
        y = F.gelu(self.fco1(y))
        y = F.gelu(self.fco2(y))
        y = F.gelu(self.fco3(y))
        out = self.fco4(y)  # (B, dim_out)
        return out

# ======================================================================================

model = TimeSeriesTransformer(dim_in=num_feat, dim_out=num_feat, pos_len=seq_len, group_size=groupsize, num_groups=numgroups)
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
mse = np.mean((Ysim - Ytrue)**2, axis=0)
print("Mean MSE:", np.mean(mse))
# show a bar plot of mse per feature
plt.figure()
plt.bar(range(len(mse)), mse)
plt.xlabel("Feature index")
plt.xticks(ticks=range(0, len(mse), groupsize), labels=feature_names[::groupsize], rotation=90, ha='right')
plt.ylabel("MSE")
plt.title("Mean Squared Error per Feature")
plt.grid(axis='y')
plt.show()

# a few example channels 
# get index of channels sorted by mse
sorted_indices = np.argsort(mse)[::-1]  # descending order
# get first, middle, and end of sorted_indices
example_indices = [sorted_indices[0], sorted_indices[len(sorted_indices)//2], sorted_indices[-1]]
# Create a single figure with vertically stacked subplots
fig, axes = plt.subplots(len(example_indices), 1, sharex=True, figsize=(10, 8))
for ax, idx in zip(axes, example_indices):
    ax.plot(Ytrue[:, idx], label="True")
    ax.plot(Ysim[:, idx], label="Simulated", alpha=0.7)
    ax.set_ylabel("Feature Value")
    ax.set_title(f"Feature: {feature_names[idx]} (MSE: {mse[idx]:.4f})")
    ax.legend()
    ax.grid(axis='both')
axes[-1].set_xlabel("Sample")  # Set x-label only on the last subplot
plt.tight_layout()
plt.show()