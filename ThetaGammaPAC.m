%% Step 1: Calculate Theta and Gamma Signals
% Run the MemoryThetaMeasurements script to calculate theta signals
run('MemoryThetaMeasurements.m'); % Ensure this script outputs theta signals
thetaSignals = {}; % Initialize cell array for theta signals

% Extract theta signals from the results structure
for idx = 1:length(thetaPowerResults)
    % full-band instantaneous phase (for whole recording)
    thetaSignals{idx} = thetaPowerResults(idx).filteredPhase; % takes dataOneChannelSelPhase
    % per-state per-trial phase cell arrays (encoding / decoding)
    theta_encodingsignal{idx} = thetaPowerResults(idx).encodingPhaseSignal; % takes encodePhase
    theta_decodingsignal{idx} = thetaPowerResults(idx).decodingPhaseSignal; % takes decodePhase
    theta_waitingsignal{idx} = thetaPowerResults(idx).waitingPhaseSignal; % takes waitPhase
    % expected - bar plots would have less error 
end

% Run the MemoryGammaMeasurements script to calculate gamma signals
run('MemoryGammaMeasurements.m'); % Ensure this script outputs gamma signals
gammaSignals = {}; % Initialize cell array for gamma signals

% Extract gamma signals from the results structure
for idx = 1:length(gammaPowerResults)
    % compute full-trace gamma amplitude envelope from filteredSignal
    %gammaSignals{idx} = abs(hilbert(gammaPowerResults(idx).filteredSignal)); - COMMENTED OUT FOR NOW
    % per-state per-trial amplitude cell arrays (encoding / decoding)
    gamma_encodingsignal{idx} = gammaPowerResults(idx).encodingAmpSignalLower; % takes encodeAmp_lower TAKES LOWER ENCODING!!!!
    gamma_decodingsignal{idx} = gammaPowerResults(idx).decodingAmpSignalLower; % takes decodeAmp_lower TAKES LOWER DECODING!!!!
    gamma_waitingsignalLower{idx} = gammaPowerResults(idx).waitingAmpSignalLower; % takes waitAmp_lower
    gamma_waitingsignalHigher{idx} = gammaPowerResults(idx).waitingAmpSignalHigher; % takes waitAmp_higher
    % PAC_encoding = zeros(size(gamma_encodingsignal));
    % for each trial split the 20x1 cell array into 20 arrays, run CalcPac
    % for each - 20 Calcpacs/channel/enc-dec. For each channel store
    % average and stdev and make a bar plot at the end based on avg val for
    % each channel
    %{
    for trl = 1:length(PAC_encoding)
        PAC_encoding(trl) = CalcPac % pass gamma and theta cell here
    %}
end

%% Step 2: Calculate PAC for Each Channel - using filtered signal
numChannels = min(length(thetaSignals), length(gammaSignals)); % Ensure matching number of channels
modulationIndices = zeros(1, numChannels); % To store MI for each channel
phaseAmplitudePlots = cell(1, numChannels); % To store plots if needed

for ch = 1:numChannels
    % Get the theta and gamma signals for the current channel
    thetaSignal = thetaSignals{ch}; % Theta phase signal
    gammaSignal = gammaSignals{ch}; % Gamma amplitude signal
    
    % Calculate PAC using calcPAChelper (expects phase, amplitude)
    [MI, P, bcent, fig1] = calcPAChelper(thetaSignal, gammaSignal, 18, true); % Enable plotting
    
    % Store the Modulation Index (MI)
    modulationIndices(ch) = MI;
    
    % Optionally, store the figure handle if you want to save or analyze later
    phaseAmplitudePlots{ch} = fig1;
    
    % Display the result for the current channel
    disp(['Channel ', num2str(ch), ' - Modulation Index (MI): ', num2str(MI)]);
end

%% Step 3: Display All Modulation Indices
disp('Modulation Indices for all channels:');
disp(modulationIndices);

% Plot all MIs in a bar chart
figure;
bar(modulationIndices, 'FaceColor', [0.2, 0.6, 0.8]);
xlabel('Channel');
ylabel('Modulation Index (MI)');
title('Phase-Amplitude Coupling (PAC) Across Channels');

%% Step 4: Calculate PAC for Encoding and Decoding Periods Separately
numChannels = min(length(theta_encodingsignal), length(gamma_encodingsignal)); 
numTrials = 20; % Number of trials per channel

% Initialize storage for PAC values
PAC_enc_all = cell(1, numChannels); % Store PAC values for all trials per channel (encoding)
PAC_dec_all = cell(1, numChannels); % Store PAC values for all trials per channel (decoding)
PAC_wait_all = cell(1, numTrials); % Store PAC values for all trials per channel (waiting) GOING TO USE LOWER HERE

% Initialize storage for average and stdev
modulationIndices_enc_avg = zeros(1, numChannels);
modulationIndices_enc_std = zeros(1, numChannels);
modulationIndices_dec_avg = zeros(1, numChannels);
modulationIndices_dec_std = zeros(1, numChannels);
modulationIndices_wait_avg = zeros(1, numChannels);
modulationIndices_wait_std = zeros(1, numChannels);

for ch = 1:numChannels
    % Get the encoding and decoding signals for the current channel (20x1 cell arrays)
    thetaSignal_enc_trials = theta_encodingsignal{ch}; % 20x1 cell array
    gammaSignal_enc_trials = gamma_encodingsignal{ch}; % 20x1 cell array
    thetaSignal_dec_trials = theta_decodingsignal{ch}; % 20x1 cell array
    gammaSignal_dec_trials = gamma_decodingsignal{ch}; % 20x1 cell array
    thetaSignal_wait_trials = theta_waitingsignal{ch}; % 20x1 cell array
    gammaSignal_wait_trials = gamma_waitingsignalLower{ch}; % 20x1 cell array - USING LOWER BAND HERE
    
    % Initialize storage for PAC values for this channel
    PAC_enc_trials = zeros(1, numTrials);
    PAC_dec_trials = zeros(1, numTrials);
    PAC_wait_trials = zeros(1, numTrials);
    
    % Loop through each trial
    for trl = 1:numTrials
        % Extract signals for the current trial
        thetaSignal_enc = thetaSignal_enc_trials{trl};
        gammaSignal_enc = gammaSignal_enc_trials{trl};
        thetaSignal_dec = thetaSignal_dec_trials{trl};
        gammaSignal_dec = gammaSignal_dec_trials{trl};
        thetaSignal_wait = thetaSignal_wait_trials{trl};
        gammaSignal_wait = gammaSignal_wait_trials{trl};
        
        % Skip if any signal is empty
        if isempty(thetaSignal_enc) || isempty(gammaSignal_enc)
            warning(['Channel ', num2str(ch), ', Trial ', num2str(trl), ' - Encoding signals are empty. Skipping...']);
            PAC_enc_trials(trl) = NaN; % Mark as NaN
            continue;
        end
        if isempty(thetaSignal_dec) || isempty(gammaSignal_dec)
            warning(['Channel ', num2str(ch), ', Trial ', num2str(trl), ' - Decoding signals are empty. Skipping...']);
            PAC_dec_trials(trl) = NaN; % Mark as NaN
            continue;
        end
        if isempty(thetaSignal_wait) || isempty(gammaSignal_wait)
            warning(['Channel ', num2str(ch), ', Trial ', num2str(trl), ' - Waiting signals are empty. Skipping...']);
            PAC_wait_trials(trl) = NaN; % Mark as NaN
            continue;
        end
        
        % Calculate PAC for Encoding using theta phase & gamma amplitude
        [MI_enc, ~, ~, ~] = calcPAChelper(thetaSignal_enc, gammaSignal_enc, 18, false); % No plotting
        PAC_enc_trials(trl) = MI_enc;
        disp(['Channel ', num2str(ch), ', Trial ', num2str(trl), ' - Encoding MI: ', num2str(MI_enc)]);
        
        % Calculate PAC for Decoding
        [MI_dec, ~, ~, ~] = calcPAChelper(thetaSignal_dec, gammaSignal_dec, 18, false); % No plotting
        PAC_dec_trials(trl) = MI_dec;
        disp(['Channel ', num2str(ch), ', Trial ', num2str(trl), ' - Decoding MI: ', num2str(MI_dec)]);

        % Calculate PAC for Waiting
        [MI_wait, ~, ~, ~] = calcPAChelper(thetaSignal_wait, gammaSignal_wait, 18, false); % No plotting
        PAC_wait_trials(trl) = MI_wait;
        disp(['Channel ', num2str(ch), ', Trial ', num2str(trl), ' - Waiting MI: ', num2str(MI_wait)]);
    end
    
    % Store PAC values for this channel
    PAC_enc_all{ch} = PAC_enc_trials;
    PAC_dec_all{ch} = PAC_dec_trials;
    PAC_wait_all{ch} = PAC_wait_trials;
    
    % Calculate average and standard deviation (ignoring NaN values)
    modulationIndices_enc_avg(ch) = mean(PAC_enc_trials, 'omitnan');
    modulationIndices_enc_std(ch) = std(PAC_enc_trials, 'omitnan');
    modulationIndices_dec_avg(ch) = mean(PAC_dec_trials, 'omitnan');
    modulationIndices_dec_std(ch) = std(PAC_dec_trials, 'omitnan');
    modulationIndices_wait_avg(ch) = mean(PAC_wait_trials, 'omitnan');
    modulationIndices_wait_std(ch) = std(PAC_wait_trials, 'omitnan');
    
    % Display the results for the current channel
    disp(['Channel ', num2str(ch), ' - Encoding MI (avg ± std): ', ...
          num2str(modulationIndices_enc_avg(ch)), ' ± ', num2str(modulationIndices_enc_std(ch))]);
    disp(['Channel ', num2str(ch), ' - Decoding MI (avg ± std): ', ...
          num2str(modulationIndices_dec_avg(ch)), ' ± ', num2str(modulationIndices_dec_std(ch))]);
    disp(['Channel ', num2str(ch), ' - Waiting MI (avg ± std): ', ...
          num2str(modulationIndices_wait_avg(ch)), ' ± ', num2str(modulationIndices_wait_std(ch))]);
end

% Validation: Display number of valid measurements per channel
disp('--- Validation: Number of Valid Measurements per Channel ---');
for ch = 1:numChannels
    numValid_enc = sum(~isnan(PAC_enc_all{ch})); % Count non-NaN encoding measurements
    numValid_dec = sum(~isnan(PAC_dec_all{ch})); % Count non-NaN decoding measurements
    numValid_wait = sum(~isnan(PAC_wait_all{ch})); % Count non-NaN waiting measurements

    disp(['Channel ', num2str(ch), ':']);
    disp(['  - Encoding: ', num2str(numValid_enc), ' / ', num2str(numTrials), ' valid measurements']);
    disp(['  - Decoding: ', num2str(numValid_dec), ' / ', num2str(numTrials), ' valid measurements']);
    disp(['  - Waiting: ', num2str(numValid_wait), ' / ', num2str(numTrials), ' valid measurements']);
end
disp('--- End of Validation ---');
%% Step 5: Display and Plot Encoding and Decoding Modulation Indices
disp('Encoding Modulation Indices (avg) for all channels:');
disp(modulationIndices_enc_avg);
disp('Decoding Modulation Indices (avg) for all channels:');
disp(modulationIndices_dec_avg);
disp('Waiting Modulation Indices (avg) for all channels:');
disp(modulationIndices_wait_avg);

% Plot Encoding and Decoding MIs in a grouped bar chart with error bars
figure;
barVals = [modulationIndices_enc_avg; modulationIndices_dec_avg; modulationIndices_wait_avg]'; % Grouped bar values
errorVals = [modulationIndices_enc_std; modulationIndices_dec_std; modulationIndices_wait_std]'; % Error bars

% Create grouped bar chart
b = bar(barVals, 'grouped');
hold on;

% Add error bars
numGroups = size(barVals, 1);
numBars = size(barVals, 2);
groupWidth = min(0.8, numBars/(numBars + 1.5));

for i = 1:numBars
    x = (1:numGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*numBars);
    errorbar(x, barVals(:,i), errorVals(:,i), 'k', 'linestyle', 'none');
end

hold off;

% Set axis properties
set(gca, 'XTick', 1:numChannels);
set(gca, 'XTickLabel', arrayfun(@num2str, 1:numChannels, 'UniformOutput', false));
legend({'Encoding', 'Decoding', 'Waiting'});
xlabel('Channel');
ylabel('Modulation Index (MI)');
title('PAC Modulation Indices for Encoding vs. Decoding vs. Waiting (Avg ± Std)');

%% Step 6: Bar Plot of Channel 1 and 2 Trial Results
channelsToPlot = [1, 2]; % Channels to plot
numTrials = 20; % Number of trials
figure;
for i = 1:length(channelsToPlot)
    ch = channelsToPlot(i);
    subplot(1, length(channelsToPlot), i);
    bar([PAC_enc_all{ch}; PAC_dec_all{ch}; PAC_wait_all{ch}]', 'grouped');
    set(gca, 'XTick', 1:numTrials);
    xlabel('Trial');
    ylabel('Modulation Index (MI)');
    title(['Channel ', num2str(ch), ' - Encoding vs Decoding MI']);
    legend({'Encoding', 'Decoding'});
end