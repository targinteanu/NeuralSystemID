import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset, Dataset
from torch.nn import GELU
import pandas as pd
from sklearn.model_selection import train_test_split

# Step 1: Define the Neural Network Architecture
class NeuralNetwork(nn.Module):
    def __init__(self, input_size, hidden_size, output_size):
        super(NeuralNetwork, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(input_size, hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, output_size),
        )

    def forward(self, x):
        return self.net(x)

# Step 2: Initialize the Model, Loss Function, and Optimizer
batch_size = 32
input_size = 101 # todo: this should be defined based width of X instead
hidden_size = 128
output_size = 100 # todo: this should be defined based on width of y instead
model = NeuralNetwork(input_size, hidden_size, output_size)

criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

# Step 3: Prepare the Data ---------------------------------------------------------------------
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

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

# Create TensorDatasets
train_dataset = TensorDataset(X_train, Y_train)
test_dataset = TensorDataset(X_test, Y_test)

# Create DataLoaders for batching
train_loader = DataLoader(train_dataset, batch_size, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size)

# -----------------------------------------------------------------------------------------------

# Step 4: Train the Model
model.train()
num_epochs = 25
for epoch in range(num_epochs):
    for X_batch, Y_batch in train_loader:
        # Forward pass
        Y_pred = model(X_batch)
        loss = criterion(Y_pred, Y_batch)

        # Backward pass and optimization
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    print(f"Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}")

# Step 5: Evaluate the Model on Test Data
model.eval()
with torch.no_grad():
    total_loss = 0
    for X_batch, Y_batch in test_loader:
        Y_pred = model(X_batch)
        loss = criterion(Y_pred, Y_batch)
        total_loss += loss.item() * X_batch.size(0)  # sum up batch loss

    avg_loss = total_loss / len(test_dataset)
    print(f"Test Loss: {avg_loss:.4f}")

# Step 6: Save the Trained Model
# torch.save(model.state_dict(), "neural_network.pth")
# print("Model saved to 'neural_network.pth'")