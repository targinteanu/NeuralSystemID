import torch
import copy
from torch.amp import GradScaler, autocast

def trainDynsysModel(
        model, 
        optimizer, 
        criterion, 
        TrainLoader, # can be single loader or tuple list to run sequentially
        TestLoader, # validation loader for early stopping
        num_epochs = 100, # max 
        patience = 10, # epochs to wait for improvement before stopping
        allow_early_stopping = True,
):
    
    # determine GPU or CPU device: 
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    model.to(device)
    
    # make a variable to store the best model
    best_model = model.state_dict()
    
    model.train()
    best_val = float('inf')
    no_improve = 0
    train_losses = []
    val_losses = []
    scaler = GradScaler()

    if not isinstance(TrainLoader, tuple):
        # TO DO: think about organizing this as dict instead of list/tuple 
        TrainLoader = (TrainLoader,)
        TestLoader = (TestLoader,)
    
    for train_loader, test_loader in zip(TrainLoader, TestLoader):

        train_size = len(train_loader.dataset)

        for epoch in range(num_epochs):
            # --- train ---
            model.train()
            running_loss = 0.0
            for X_batch, Y_batch, U_batch in train_loader:
                X_batch, Y_batch, U_batch = X_batch.to(device), Y_batch.to(device), U_batch.to(device)
                with autocast(device_type='cuda' if torch.cuda.is_available() else 'cpu'):
                    # Forward pass
                    Y_pred = model(X_batch, U_batch)
                    loss = criterion(Y_pred, Y_batch)
                # Backward pass and optimization
                optimizer.zero_grad()
                scaler.scale(loss).backward()
                scaler.step(optimizer)
                scaler.update()
                running_loss += loss.item() * X_batch.size(0)
            epoch_train_loss = running_loss / train_size
            train_losses.append(epoch_train_loss)

            # --- validate ---
            model.eval()
            val_running = 0.0
            with torch.no_grad():
                for X_val, Y_val, U_val in test_loader:   # use test_loader or a separate val_loader
                    X_val, Y_val, U_val = X_val.to(device), Y_val.to(device), U_val.to(device)
                    with autocast(device_type='cuda' if torch.cuda.is_available() else 'cpu'):
                        Y_val_pred = model(X_val, U_val)
                        l = criterion(Y_val_pred, Y_val)
                    val_running += l.item() * X_val.size(0)
            epoch_val_loss = val_running / len(test_loader.dataset)
            val_losses.append(epoch_val_loss)

            print(f"Epoch {epoch+1}/{num_epochs} — train_loss: {epoch_train_loss:.6f}, val_loss: {epoch_val_loss:.6f}")

            # --- Early stopping and best model recording ---
            if allow_early_stopping:
                if epoch_val_loss < best_val - 1e-9:
                    best_val = epoch_val_loss
                    best_model = copy.deepcopy(model.state_dict())  # save the best model weights
                    no_improve = 0
                else:
                    no_improve += 1
                    if no_improve >= patience:
                        print(f"No improvement for {patience} epochs — stopping early at epoch {epoch+1}.")
                        break

        # Load the best model weights before returning
        if allow_early_stopping:
            print(f"Restoring the best model with val_loss: {best_val:.6f}")
            missing, unexpected = model.load_state_dict(best_model)
            if missing:
                print(f"Warning: Missing keys when loading best model: {missing}")
            if unexpected:
                print(f"Warning: Unexpected keys when loading best model: {unexpected}")
        else:
            print("Early stopping disabled. Returning the final model.")

    return model, train_losses, val_losses