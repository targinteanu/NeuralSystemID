import pandas as pd
import numpy as np
import torch
from torch.utils.data import Dataset

# ------------------------
# Load main data
# ------------------------
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
import bisect

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
event_counts = np.array([count_events_in_window(t) for t in time])

# Append event_counts as an additional input feature
data_aug = np.hstack([data, event_counts.reshape(-1, 1)])
# Now each input row has: [original data..., event_count]

# ------------------------
# Build input-output pairs using the 3-row / 3 ms rule
# ------------------------
dt_target = 0.003      # 3 ms
dt_tol = 0.0005        # ±0.5 ms

inputs = []
outputs = []

N = len(df)

for i in range(N - 3):
    dt = time[i+3] - time[i]
    if abs(dt - dt_target) <= dt_tol:
        inputs.append(data_aug[i])   # augmented input with event_count
        outputs.append(data[i+3])    # output is ONLY the data (no event count)

X = torch.tensor(inputs, dtype=torch.float32)
Y = torch.tensor(outputs, dtype=torch.float32)

print("Pairs created:", len(X))
print("Input shape :", X.shape)   # features + 1
print("Output shape:", Y.shape)

# ------------------------
# Optional PyTorch Dataset
# ------------------------
class TimeShiftDataset(Dataset):
    def __init__(self, X, Y):
        self.X = X
        self.Y = Y

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.Y[idx]

dataset = TimeShiftDataset(X, Y)