import pandas as pd
import numpy as np
import math
import bisect


def prepTimeData(
    dt_target = 0.005, # model sample time, s
    eventwindow = 0.2, # sliding time window for counting events, s
    datafile = "data.csv",
    baselinedatafile = "baselinedata.csv",
    eventsfile = "events.csv", 
    infofile = "info.csv", 
    datarawfile = "dataraw.csv",
    baselinedatarawfile = "baselinedataraw.csv",
    maxNumel = 2**30, # max number of elements in X
    drawFromStart = True, # pull samples from beginning vs end
    omitOutliers = True,
):
    """
    Prepare input-output data from CSV files for model training.
    Loads input data "X" into shape (num_pairs, num_features)
    and output data "Y" into shape (num_pairs, num_outputs).
    Should be equivalent to prepTimeSeqData with seq_len=1.
    """

    # ------------------------
    # Load info
    # ------------------------
    info_df = pd.read_csv(infofile)
    fs = info_df.iloc[0, 5]                  # sampling frequency (Hz)
    print("Sampling frequency (Hz):", fs)
    feature_names = info_df.iloc[:, 0].values
    maxPairs = math.floor(maxNumel / (len(feature_names)))

    # ------------------------
    # Load raw data
    # ------------------------
    baseline_df = pd.read_csv(baselinedatarawfile)
    baseline_time = baseline_df.iloc[:, 0].values            # time column (seconds)
    baseline_data_raw = baseline_df.iloc[:, 1:].values           # data columns

    if drawFromStart:
        iBLend = min(baseline_data_raw.shape[0], maxPairs)
        baseline_time = baseline_time[:iBLend]
        baseline_data_raw = baseline_data_raw[:iBLend, :]
    else: 
        iBLstart = max(0, baseline_data_raw.shape[0] - maxPairs)
        baseline_time = baseline_time[iBLstart:]
        baseline_data_raw = baseline_data_raw[iBLstart:, :]

    df = pd.read_csv(datarawfile)
    time = df.iloc[:, 0].values            # time column (seconds)
    data_raw = df.iloc[:, 1:].values           # data columns

    if drawFromStart:
        iEnd = min(data_raw.shape[0], maxPairs)
        time = time[:iEnd]
        data_raw = data_raw[:iEnd, :]
    else: 
        iStart = max(0, data_raw.shape[0] - maxPairs)
        time = time[iStart:]
        data_raw = data_raw[iStart:, :]

    # ------------------------
    # Determine outliers
    # ------------------------
    if omitOutliers:
        threshsd = 3 # standard deviations 
        threshprop = .5 # proportion of features
        BLmean = np.mean(baseline_data_raw, axis=0)
        BLstd = np.std(baseline_data_raw, axis=0)
        BLisout = np.abs(baseline_data_raw - BLmean) > (threshsd * BLstd)
        BLisnoise = np.sum(BLisout, axis=1) > (threshprop * baseline_data_raw.shape[1])
        isout = np.abs(data_raw - BLmean) > (threshsd * BLstd)
        isnoise = np.sum(isout, axis=1) > (threshprop * data_raw.shape[1])
    else:
        BLisnoise = np.zeros(baseline_data_raw.shape[0], dtype=bool)
        isnoise = np.zeros(data_raw.shape[0], dtype=bool)

    # ------------------------
    # Load processed data
    # ------------------------
    baseline_df = pd.read_csv(baselinedatafile)
    baseline_time = baseline_df.iloc[:, 0].values            # time column (seconds)
    baseline_data = baseline_df.iloc[:, 1:].values           # data columns

    if drawFromStart:
        baseline_time = baseline_time[:iBLend]
        baseline_data = baseline_data[:iBLend, :]
    else:
        baseline_time = baseline_time[iBLstart:]
        baseline_data = baseline_data[iBLstart:, :]

    df = pd.read_csv(datafile)
    time = df.iloc[:, 0].values            # time column (seconds)
    data = df.iloc[:, 1:].values           # data columns

    if drawFromStart:
        time = time[:iEnd]
        data = data[:iEnd, :]
    else:
        time = time[iStart:]
        data = data[iStart:, :]

    # ------------------------
    # Load events
    # ------------------------
    events_df = pd.read_csv(eventsfile)
    event_times = events_df.iloc[:, 0].values  # assume first column is event time in seconds
    event_times = np.sort(event_times)         # ensure sorted

    # ------------------------
    # Compute event count for each row
    # ------------------------
    # For fast lookup using binary search
    def count_events_in_window(t, window=0.2):
        """
        Count how many event_times fall in (t - window, t].
        Uses bisect for O(log n) search.
        """
        left = bisect.bisect_right(event_times, t - window)
        right = bisect.bisect_right(event_times, t)
        return right - left

    event_counts = np.array([count_events_in_window(t, eventwindow) for t in time])
    event_counts = event_counts.reshape(-1, 1)

    # Append event_counts as an additional input feature
    data_aug = np.hstack([data, event_counts])
    data_aug_raw = np.hstack([data_raw, event_counts])
    # Now each input row has: [original data..., event_count]
    baseline_data_aug = np.hstack([baseline_data, np.zeros((baseline_data.shape[0], 1))]) 
    baseline_data_aug_raw = np.hstack([baseline_data_raw, np.zeros((baseline_data_raw.shape[0], 1))])
    # baseline has zero event counts

    # ------------------------
    # Build input-output pairs using the _ ms rule
    # ------------------------
    dt_tol = 0.15 * dt_target
    drow_target = int(dt_target * fs)  # number of rows 

    # ------------------------
    # baseline data 
    # ------------------------
    X_list = []
    Y_list = []
    Xr_list = []
    Yr_list = []

    if drawFromStart:
        Nbl = iBLend
    else:
        Nbl = len(baseline_data_aug)
    for i in range(Nbl - drow_target):
        dt = baseline_time[i+drow_target] - baseline_time[i]
        if (abs(dt - dt_target) <= dt_tol) and (not BLisnoise[i]) and (not BLisnoise[i+drow_target]):
            X_list.append(baseline_data_aug[i])   
            Y_list.append(baseline_data[i+drow_target])    
            Xr_list.append(baseline_data_aug_raw[i])
            Yr_list.append(baseline_data_raw[i+drow_target])

    Xbl = np.concatenate(X_list, axis=0)
    Ybl = np.concatenate(Y_list, axis=0)
    Xbl_raw = np.concatenate(Xr_list, axis=0)
    Ybl_raw = np.concatenate(Yr_list, axis=0)

    print("Baseline Pairs created:", len(Xbl))
    print("Input shape :", Xbl.shape)   
    print("Output shape:", Ybl.shape)

    # ------------------------
    # main data 
    # ------------------------
    X_list = []
    Y_list = []
    Xr_list = []
    Yr_list = []

    if drawFromStart:
        N = iEnd
    else:
        N = len(data_aug)
    for i in range(N - drow_target):
        dt = time[i+drow_target] - time[i]
        if (abs(dt - dt_target) <= dt_tol) and (not isnoise[i]) and (not isnoise[i+drow_target]):
            X_list.append(data_aug[i])   
            Y_list.append(data[i+drow_target])    
            Xr_list.append(data_aug_raw[i])
            Yr_list.append(data_raw[i+drow_target])

    X = np.concatenate(X_list, axis=0)
    Y = np.concatenate(Y_list, axis=0)
    X_raw = np.concatenate(Xr_list, axis=0)
    Y_raw = np.concatenate(Yr_list, axis=0)

    print("Main Pairs created:", len(X))
    print("Input shape :", X.shape)   
    print("Output shape:", Y.shape)

    return fs, feature_names, X, Y, Xbl, Ybl, X_raw, Y_raw, Xbl_raw, Ybl_raw


def prepTimeSeqData(
    dt_target = 0.005, # model sample time, s
    seq_len = 64, # model transformer samples
    hzn_len = 1, # samples
    datafile = "data.csv",
    baselinedatafile = "baselinedata.csv",
    eventsfile = "events.csv", 
    infofile = "info.csv", 
    datarawfile = "dataraw.csv",
    baselinedatarawfile = "baselinedataraw.csv",
    maxNumel = 2**30, # max number of elements in X
    drawFromStart = True, # pull samples from beginning vs end
    omitOutliers = True,
):
    """
    Prepare time-sequence data from CSV files for model training.
    Loads input data "X" into shape (num_pairs, seq_len, num_features)
    and output data "Y" into shape (num_pairs, num_outputs).
    """

    # ------------------------
    # Load info
    # ------------------------
    info_df = pd.read_csv(infofile)
    fs = info_df.iloc[0, 5]                  # sampling frequency (Hz)
    print("Sampling frequency (Hz):", fs)
    feature_names = info_df.iloc[:, 0].values
    maxPairs = math.floor(maxNumel / (seq_len * len(feature_names)))

    # ------------------------
    # Load raw data
    # ------------------------
    baseline_df = pd.read_csv(baselinedatarawfile)
    baseline_time = baseline_df.iloc[:, 0].values            # time column (seconds)
    baseline_data_raw = baseline_df.iloc[:, 1:].values           # data columns

    if drawFromStart:
        iBLend = min(baseline_data_raw.shape[0], maxPairs)
        baseline_time = baseline_time[:iBLend]
        baseline_data_raw = baseline_data_raw[:iBLend, :]
    else: 
        iBLstart = max(0, baseline_data_raw.shape[0] - maxPairs)
        baseline_time = baseline_time[iBLstart:]
        baseline_data_raw = baseline_data_raw[iBLstart:, :]

    df = pd.read_csv(datarawfile)
    time = df.iloc[:, 0].values            # time column (seconds)
    data_raw = df.iloc[:, 1:].values           # data columns

    if drawFromStart:
        iEnd = min(data_raw.shape[0], maxPairs)
        time = time[:iEnd]
        data_raw = data_raw[:iEnd, :]
    else: 
        iStart = max(0, data_raw.shape[0] - maxPairs)
        time = time[iStart:]
        data_raw = data_raw[iStart:, :]

    # ------------------------
    # Determine outliers
    # ------------------------
    if omitOutliers:
        threshsd = 3 # standard deviations 
        threshprop = .5 # proportion of features
        BLmean = np.mean(baseline_data_raw, axis=0)
        BLstd = np.std(baseline_data_raw, axis=0)
        BLisout = np.abs(baseline_data_raw - BLmean) > (threshsd * BLstd)
        BLisnoise = np.sum(BLisout, axis=1) > (threshprop * baseline_data_raw.shape[1])
        isout = np.abs(data_raw - BLmean) > (threshsd * BLstd)
        isnoise = np.sum(isout, axis=1) > (threshprop * data_raw.shape[1])
    else:
        BLisnoise = np.zeros(baseline_data_raw.shape[0], dtype=bool)
        isnoise = np.zeros(data_raw.shape[0], dtype=bool)

    # ------------------------
    # Load processed data
    # ------------------------
    baseline_df = pd.read_csv(baselinedatafile)
    baseline_time = baseline_df.iloc[:, 0].values            # time column (seconds)
    baseline_data = baseline_df.iloc[:, 1:].values           # data columns

    if drawFromStart:
        baseline_time = baseline_time[:iBLend]
        baseline_data = baseline_data[:iBLend, :]
    else:
        baseline_time = baseline_time[iBLstart:]
        baseline_data = baseline_data[iBLstart:, :]

    df = pd.read_csv(datafile)
    time = df.iloc[:, 0].values            # time column (seconds)
    data = df.iloc[:, 1:].values           # data columns

    if drawFromStart:
        time = time[:iEnd]
        data = data[:iEnd, :]
    else:
        time = time[iStart:]
        data = data[iStart:, :]

    # ------------------------
    # Load events
    # ------------------------
    events_df = pd.read_csv(eventsfile)
    event_times = events_df.iloc[:, 0].values  # assume first column is event time in seconds
    event_times = np.sort(event_times)         # ensure sorted

    # ------------------------
    # Compute event count for each row
    # ------------------------
    # For fast lookup using binary search
    def count_events_in_window(t, window=0.2):
        """
        Count how many event_times fall in (t - window, t].
        Uses bisect for O(log n) search.
        """
        left = bisect.bisect_right(event_times, t - window)
        right = bisect.bisect_right(event_times, t)
        return right - left

    event_counts = np.array([count_events_in_window(t, 1.2/fs) for t in time])
    event_counts = event_counts.reshape(-1, 1)

    # Append event_counts as an additional input feature
    data_aug = np.hstack([data, event_counts])
    data_aug_raw = np.hstack([data_raw, event_counts])
    # Now each input row has: [original data..., event_count]
    baseline_data_aug = np.hstack([baseline_data, np.zeros((baseline_data.shape[0], 1))]) 
    baseline_data_aug_raw = np.hstack([baseline_data_raw, np.zeros((baseline_data_raw.shape[0], 1))]) 
    # baseline has zero event counts

    # ------------------------
    # Build input-output pairs using the _ ms rule
    # ------------------------
    dt_tol = 0.15 * dt_target
    drow_target = int(dt_target * fs)  # number of rows 

    # create sliding windows
    def create_windows(data, data2, seq_len=128, horizon=1):
        X, Y, W, Z = [], [], [], []
        for i in range(len(data) - seq_len - horizon + 1):
            X.append(data[i:i+seq_len])
            Y.append(data[i+seq_len + horizon - 1])
            W.append(data2[i:i+seq_len])
            Z.append(data2[i+seq_len + horizon - 1])
        return np.array(X), np.array(Y), np.array(W), np.array(Z)

    # ------------------------
    # baseline data 
    # ------------------------
    X_list = []
    Y_list = []
    Xr_list = []
    Yr_list = []

    if drawFromStart:
        Nbl = iBLend
    else:
        Nbl = len(baseline_data_aug)
    for iStart in range(drow_target):
        inputs = []
        inputs_raw = []
        for i in range(iStart, (Nbl - drow_target), drow_target):
            dt = baseline_time[i+drow_target] - baseline_time[i]
            if (abs(dt - dt_target) <= dt_tol) and (not BLisnoise[i]):
                inputs.append(baseline_data_aug[i, :]) 
                inputs_raw.append(baseline_data_aug_raw[i, :])
            else:
                if len(inputs) > seq_len+hzn_len:
                    x, y, xr, yr = create_windows(inputs, inputs_raw, seq_len, hzn_len)
                    if x is not None:
                        X_list.append(x)
                        Y_list.append(y)
                        Xr_list.append(xr)
                        Yr_list.append(yr)
                inputs = []
                inputs_raw = []
        # catch trailing segment
        if len(inputs) > seq_len+hzn_len:
            x, y, xr, yr = create_windows(inputs, inputs_raw, seq_len, hzn_len)
            if x is not None:
                X_list.append(x)
                Y_list.append(y)
                Xr_list.append(xr)
                Yr_list.append(yr)

    Xbl = np.concatenate(X_list, axis=0)
    Ybl = np.concatenate(Y_list, axis=0)
    Xbl_raw = np.concatenate(Xr_list, axis=0)
    Ybl_raw = np.concatenate(Yr_list, axis=0)

    print("Baseline Pairs created:", len(Xbl))
    print("Input shape :", Xbl.shape)   
    print("Output shape:", Ybl.shape)

    # ------------------------
    # main data 
    # ------------------------
    X_list = []
    Y_list = []
    Xr_list = []
    Yr_list = []

    if drawFromStart:
        N = iEnd
    else:
        N = len(data_aug)
    for iStart in range(drow_target):
        inputs = []
        inputs_raw = []
        for i in range(iStart, (N - drow_target), drow_target):
            dt = time[i+drow_target] - time[i]
            if (abs(dt - dt_target) <= dt_tol) and (not BLisnoise[i]):
                inputs.append(data_aug[i, :]) 
                inputs_raw.append(data_aug_raw[i, :])
            else:
                if len(inputs) > seq_len+hzn_len:
                    x, y, xr, yr = create_windows(inputs, inputs_raw, seq_len, hzn_len)
                    if x is not None:
                        X_list.append(x)
                        Y_list.append(y)
                        Xr_list.append(xr)
                        Yr_list.append(yr)
                inputs = []
                inputs_raw = []
        # catch trailing segment
        if len(inputs) > seq_len+hzn_len:
            x, y, xr, yr = create_windows(inputs, inputs_raw, seq_len, hzn_len)
            if x is not None:
                X_list.append(x)
                Y_list.append(y)
                Xr_list.append(xr)
                Yr_list.append(yr)

    X = np.concatenate(X_list, axis=0)
    Y = np.concatenate(Y_list, axis=0)
    X_raw = np.concatenate(Xr_list, axis=0)
    Y_raw = np.concatenate(Yr_list, axis=0)

    print("Main Pairs created:", len(X))
    print("Input shape :", X.shape)   
    print("Output shape:", Y.shape)

    return fs, feature_names, X, Y, Xbl, Ybl, X_raw, Y_raw, Xbl_raw, Ybl_raw