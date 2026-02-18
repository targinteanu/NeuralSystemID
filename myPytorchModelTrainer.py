import torch

def trainDynsysModel(
        model, 
        optimizer, 
        criterion, 
        train_loader,
        test_loader, # validation loader for early stopping
        num_epochs = 100, # max 
        patience = 10, # epochs to wait for improvement before stopping
        allow_early_stopping = True,
):
    
    # make a variable to store the best model
    best_model = model.state_dict()
    
    model.train()
    best_val = float('inf')
    no_improve = 0
    train_losses = []
    val_losses = []
    train_size = len(train_loader.dataset)

    for epoch in range(num_epochs):
        # --- train ---
        model.train()
        running_loss = 0.0
        for X_batch, Y_batch, U_batch in train_loader:
            # Forward pass
            Y_pred = model(X_batch, U_batch)
            loss = criterion(Y_pred, Y_batch)
            # Backward pass and optimization
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            running_loss += loss.item() * X_batch.size(0)
        epoch_train_loss = running_loss / train_size
        train_losses.append(epoch_train_loss)

        # --- validate ---
        model.eval()
        val_running = 0.0
        with torch.no_grad():
            for X_val, Y_val, U_val in test_loader:   # use test_loader or a separate val_loader
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
                best_model = model.state_dict()  # save the best model weights
                no_improve = 0
            else:
                no_improve += 1
                if no_improve >= patience:
                    print(f"No improvement for {patience} epochs — stopping early at epoch {epoch+1}.")
                    break

    # Load the best model weights before returning
    if allow_early_stopping:
        model.load_state_dict(best_model)
        
    return model, train_losses, val_losses