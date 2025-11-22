%% Calculate and display Gamma measurements 
% Adapted from MemoryThetaMeasurements.m

gammaPowerResults = struct();

%% Loading the data 
filepath = '~/Documents/Anderson Lab/Saved Data 2025-05-22 11.44.08'; 
makefullfile = @(f) fullfile(f.folder, f.name);
OnlineFiles = dir([filepath,filesep,'OnlineDisplaySavedData*.mat']);
OnlineFile = makefullfile(OnlineFiles(1)); 
NS2Files = dir([filepath,filesep,'*.ns2']); 
NS2File = makefullfile(NS2Files(1)); 
NS5Files = dir([filepath,filesep,'*.ns5']); 
NS5File = makefullfile(NS5Files(1)); 
NEVFiles = dir([filepath,filesep,'*.nev']); 
NEVFile = makefullfile(NEVFiles(1));
load(OnlineFile, "SerialLog"); 
ns2 = openNSx(NS2File, 'uV'); NS2 = ns2timetable(ns2);
ns5 = openNSx(NS5File, 'uV'); NS5 = ns2timetable(ns5);
try
    nev = openNEV(NEVFile, 'nosave', 'nomat', 'uV');
    NEV = nev2table(nev);
catch ME
    warning(ME.identifier, 'Error loading NEV file: %s', ME.message);
    NEV = []; % Set NEV to an empty array if loading fails
end
codesToKeep = [61 10 11 1:3 5 6 70 71 7 8 9]; % Default mapping codes
NEVevents = parseNEVSerialCodes(NEV, codesToKeep); % returns table with EventTime, EventCode, EventName

%% Looping through the channels 
SamplingFreq = ns2.MetaTags.SamplingFreq;
t = NS2.Time';
tRel = seconds(t-t(1));

% Select channels of interest: 
channelNames = NS2.Properties.VariableNames;
channelIndices = contains(channelNames, 'LH'); % Left Hippocampus
channelIndices = find(channelIndices);

% Stimulation times (pulse train): 
StimTrainRec = NS5.AINP1; % Analog Input 1
StimTrainRec = movmean(StimTrainRec, round(ns5.MetaTags.SamplingFreq/SamplingFreq));
StimTrainRec = interp1(NS5.Time, StimTrainRec, NS2.Time, "nearest");
StimTrainRec = [false; diff(StimTrainRec) > 2000]';

% Setup gamma bandpass filters
% Encoding range: 30-50 Hz (lower gamma)
% Decoding range: 65-100 Hz (higher gamma)
gammaBPF_lower = buildFIRBPF(SamplingFreq, 30, 50); % Encoding gamma takes from here
gammaBPF_higher = buildFIRBPF(SamplingFreq, 65, 100); % Decoding gamma takes from here

for idx = 1:length(channelIndices)

channelIndex = channelIndices(idx);
channelName = NS2.Properties.VariableNames{channelIndex}
dataOneChannel = NS2{:,channelIndex}';

dataOneChannelWithArtifact = dataOneChannel; 
t0 = t(1);

%% Artifact detection
artExtend = 10; % extend artifact by __ samples 
artIndAll = StimTrainRec; 
artIndAll = movsum(artIndAll, artExtend) > 0;
artIndAll = find(artIndAll);
[~,baselineStartInd] = max(diff(artIndAll));
baselineEndInd = artIndAll(baselineStartInd+1); baselineStartInd = artIndAll(baselineStartInd); 

%% Set baseline to fit model 
if ~isempty(artIndAll)
baselineWinLen = 1000; ARlen = 10; % samples 
dataBaseline = dataOneChannel(baselineStartInd:baselineEndInd); 
% Use the encoding-range filter for baseline AR model fitting
dataBaseline = dataBaseline - mean(dataBaseline);
baselineWin = (baselineEndInd-baselineStartInd) + [-1,1]*baselineWinLen; 
baselineWin = baselineWin/2; baselineWin = round(baselineWin); 
baselineWin(1) = max(1,baselineWin(1)); baselineWin(2) = min(length(dataBaseline),baselineWin(2));
dataBaseline = dataBaseline(baselineWin(1):baselineWin(2));
ARmdl = ar(iddata(dataBaseline', [], 1/SamplingFreq), ARlen, 'yw');

% part of artifact removal (try commenting out) - changed to remove DC offset
%% Remove artifact 
dataOneChannel = dataOneChannelWithArtifact;
dataOneChannel = dataOneChannel - mean(dataOneChannel); % correct DC offset

for ind = artIndAll
    ind0 = ind - ARlen;
    if ind0 > 0
        dataOneChannel(ind) = myFastForecastAR(ARmdl, dataOneChannel(ind0:(ind-1))', 1);
    end
end
end

% Filter raw channel twice: once for lower gamma and once for higher gamma
dataOneChannel_lower = filtfilt(gammaBPF_lower,1,dataOneChannel);
dataOneChannel_higher = filtfilt(gammaBPF_higher,1,dataOneChannel);

% analytic envelopes (instantaneous amplitude) computed over the full filtered traces
dataOneChannelAmp_lower = abs(hilbert(dataOneChannel_lower));
dataOneChannelAmp_higher = abs(hilbert(dataOneChannel_higher));

%% Select time of interest (manually)
iFirstStim = find(StimTrainRec); iFirstStim = iFirstStim(1);
selind = t <= t(iFirstStim) - seconds(60);
tSel = t(selind); tRelSel = tRel(selind);
% Create selections for each filtered trace (signal) and their amplitude envelopes
dataOneChannelSel_lower = dataOneChannel_lower(selind);
dataOneChannelSel_higher = dataOneChannel_higher(selind); % selind removes time with artifacts removed in the future we'll reconsider non-time selected

% Also select amplitude envelopes over the same time window
dataOneChannelSelAmp_lower = dataOneChannelAmp_lower(selind);
dataOneChannelSelAmp_higher = dataOneChannelAmp_higher(selind);
selind = find(selind);

%% Determine encode/decode phases of experiment 
% Use NEV serial code parser to set encode/decode/wait start/end times using default mapping

% Encoding: start = image_1, end = end_encoding
encodeStartRows = strcmpi(NEVevents.EventName, 'image_1');
encodeEndRows   = strcmpi(NEVevents.EventName, 'end_encoding');
encodeStart = NEVevents.Time(encodeStartRows);
encodeEnd   = NEVevents.Time(encodeEndRows);

% Decoding: start = decoding_start, end = decoding_end
decodeStartRows = strcmpi(NEVevents.EventName, 'decoding_start');
decodeEndRows   = strcmpi(NEVevents.EventName, 'decoding_end');
decodeStart = NEVevents.Time(decodeStartRows);
decodeEnd   = NEVevents.Time(decodeEndRows);

% Wait: start = end_encoding, end = decoding_start
waitStart = NEVevents.Time(encodeEndRows);
waitEnd   = NEVevents.Time(decodeStartRows);

%% Analyze measurements during encode/decode periods

% Extract state segments for each filter range

% Extract signal segments for each state (per-trial, per-channel)
encodeData_lower = getStateData(tSel, dataOneChannelSel_lower, encodeStart, encodeEnd);
encodeData_higher = getStateData(tSel, dataOneChannelSel_higher, encodeStart, encodeEnd);
decodeData_lower = getStateData(tSel, dataOneChannelSel_lower, decodeStart, decodeEnd);
decodeData_higher = getStateData(tSel, dataOneChannelSel_higher, decodeStart, decodeEnd);
waitData_lower = getStateData(tSel, dataOneChannelSel_lower, waitStart, waitEnd);
waitData_higher = getStateData(tSel, dataOneChannelSel_higher, waitStart, waitEnd);

% Extract amplitude-envelope segments for each state (for PAC/analysis)
encodeAmp_lower = getStateData(tSel, dataOneChannelSelAmp_lower, encodeStart, encodeEnd);
encodeAmp_higher = getStateData(tSel, dataOneChannelSelAmp_higher, encodeStart, encodeEnd);
decodeAmp_lower = getStateData(tSel, dataOneChannelSelAmp_lower, decodeStart, decodeEnd);
decodeAmp_higher = getStateData(tSel, dataOneChannelSelAmp_higher, decodeStart, decodeEnd);
waitAmp_lower = getStateData(tSel, dataOneChannelSelAmp_lower, waitStart, waitEnd);
waitAmp_higher = getStateData(tSel, dataOneChannelSelAmp_higher, waitStart, waitEnd);

% Compute Gamma Power for Encoding (30-90) and Decoding (60-120)
[avgGammaPowerEncodingLower, stdGammaPowerEncodingLower] = computeGammaPower(encodeData_lower, SamplingFreq);
[avgGammaPowerDecodingLower, stdGammaPowerDecodingLower] = computeGammaPower(decodeData_lower, SamplingFreq);
[avgGammaPowerWaitingLower, stdGammaPowerWaitingLower] = computeGammaPower(waitData_lower, SamplingFreq);
[avgGammaPowerEncodingHigher, stdGammaPowerEncodingHigher] = computeGammaPower(encodeData_higher, SamplingFreq);
[avgGammaPowerDecodingHigher, stdGammaPowerDecodingHigher] = computeGammaPower(decodeData_higher, SamplingFreq);
[avgGammaPowerWaitingHigher, stdGammaPowerWaitingHigher] = computeGammaPower(waitData_higher, SamplingFreq);

disp(['Avg Gamma Power Encoding (Lower): ', num2str(avgGammaPowerEncodingLower), ' ± ', num2str(stdGammaPowerEncodingLower)]);
disp(['Avg Gamma Power Decoding (Lower): ', num2str(avgGammaPowerDecodingLower), ' ± ', num2str(stdGammaPowerDecodingLower)]);
disp(['Avg Gamma Power Waiting (Lower): ', num2str(avgGammaPowerWaitingLower), ' ± ', num2str(stdGammaPowerWaitingLower)]);
disp(['Avg Gamma Power Encoding (Higher): ', num2str(avgGammaPowerEncodingHigher), ' ± ', num2str(stdGammaPowerEncodingHigher)]);
disp(['Avg Gamma Power Decoding (Higher): ', num2str(avgGammaPowerDecodingHigher), ' ± ', num2str(stdGammaPowerDecodingHigher)]);
disp(['Avg Gamma Power Waiting (Higher): ', num2str(avgGammaPowerWaitingHigher), ' ± ', num2str(stdGammaPowerWaitingHigher)]);

gammaPowerResults(idx).channelName = channelNames{channelIndex}; % Store channel name
gammaPowerResults(idx).encodingPowerLower = avgGammaPowerEncodingLower; % lower gamma encoding power
gammaPowerResults(idx).decodingPowerLower = avgGammaPowerDecodingLower; % lower gamma decoding power
gammaPowerResults(idx).encodingErrorLower = stdGammaPowerEncodingLower;
gammaPowerResults(idx).decodingErrorLower = stdGammaPowerDecodingLower;
gammaPowerResults(idx).encodingPowerHigher = avgGammaPowerEncodingHigher; % higher gamma encoding power
gammaPowerResults(idx).decodingPowerHigher = avgGammaPowerDecodingHigher; % higher gamma decoding power
gammaPowerResults(idx).encodingErrorHigher = stdGammaPowerEncodingHigher;
gammaPowerResults(idx).decodingErrorHigher = stdGammaPowerDecodingHigher;
gammaPowerResults(idx).waitingPowerLower = avgGammaPowerWaitingLower; % lower gamma waiting power
gammaPowerResults(idx).waitingErrorLower = stdGammaPowerWaitingLower; 
gammaPowerResults(idx).waitingPowerHigher = avgGammaPowerWaitingHigher; % higher gamma waiting power
gammaPowerResults(idx).waitingErrorHigher = stdGammaPowerWaitingHigher;

gammaPowerResults(idx).filteredSignal_lower = dataOneChannelSel_lower; % lower gamma full filtered signal
gammaPowerResults(idx).filteredSignal_higher = dataOneChannelSel_higher; % higher gamma full filtered signal

gammaPowerResults(idx).filteredAmp_lower = dataOneChannelSelAmp_lower;
gammaPowerResults(idx).filteredAmp_higher = dataOneChannelSelAmp_higher;

gammaPowerResults(idx).encodingSignalLower = encodeData_lower;
gammaPowerResults(idx).decodingSignalLower = decodeData_lower;
gammaPowerResults(idx).waitingSignalLower = waitData_lower;

gammaPowerResults(idx).encodingSignalHigher = encodeData_higher;
gammaPowerResults(idx).decodingSignalHigher = decodeData_higher;
gammaPowerResults(idx).waitingSignalHigher = waitData_higher;

gammaPowerResults(idx).encodingAmpSignalLower = encodeAmp_lower;
gammaPowerResults(idx).decodingAmpSignalLower = decodeAmp_lower;
gammaPowerResults(idx).waitingAmpSignalLower = waitAmp_lower;

gammaPowerResults(idx).encodingAmpSignalHigher = encodeAmp_higher;
gammaPowerResults(idx).decodingAmpSignalHigher = decodeAmp_higher;
gammaPowerResults(idx).waitingAmpSignalHigher = waitAmp_higher;





end

channelNamesList = {gammaPowerResults.channelName};
encodingPowersLower = [gammaPowerResults.encodingPowerLower]; 
decodingPowersLower = [gammaPowerResults.decodingPowerLower]; 
encodingErrorsLower = [gammaPowerResults.encodingErrorLower]; 
decodingErrorsLower = [gammaPowerResults.decodingErrorLower];
waitingPowersLower = [gammaPowerResults.waitingPowerLower];
waitingErrorsLower = [gammaPowerResults.waitingErrorLower];

encodingPowersHigher = [gammaPowerResults.encodingPowerHigher]; 
decodingPowersHigher = [gammaPowerResults.decodingPowerHigher]; 
encodingErrorsHigher = [gammaPowerResults.encodingErrorHigher]; 
decodingErrorsHigher = [gammaPowerResults.decodingErrorHigher];
waitingPowersHigher = [gammaPowerResults.waitingPowerHigher];
waitingErrorsHigher = [gammaPowerResults.waitingErrorHigher];

% Filter out LH03 from the results
excludeChannel = 'LH03';
includeIndices = ~strcmp(channelNamesList, excludeChannel);

% Apply the filter
channelNamesList = channelNamesList(includeIndices);
encodingPowersLower = encodingPowersLower(includeIndices);
decodingPowersLower = decodingPowersLower(includeIndices);
encodingErrorsLower = encodingErrorsLower(includeIndices);
decodingErrorsLower = decodingErrorsLower(includeIndices);
waitingPowersLower = waitingPowersLower(includeIndices);
waitingErrorsLower = waitingErrorsLower(includeIndices);

encodingPowersHigher = encodingPowersHigher(includeIndices);
decodingPowersHigher = decodingPowersHigher(includeIndices);
encodingErrorsHigher = encodingErrorsHigher(includeIndices);
decodingErrorsHigher = decodingErrorsHigher(includeIndices);
waitingPowersHigher = waitingPowersHigher(includeIndices);
waitingErrorsHigher = waitingErrorsHigher(includeIndices);
%%
figure;
barVals = [encodingPowersLower; decodingPowersLower; waitingPowersLower]'; % Grouped bar values
bar(barVals, 'grouped'); % Plot bar chart

% Compute X positions for error bars
hold on;
numGroups = size(barVals, 1);
numBars = size(barVals, 2);
groupWidth = min(0.8, numBars / (numBars + 1.5));
x = nan(numBars, numGroups);

for i = 1:numBars
    x(i,:) = (1:numGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*numBars);
end

% Add error bars
% Use the Lower-band arrays for errorbars here (previously used wrong vars)
errorbar(x(1,:), encodingPowersLower, encodingErrorsLower, 'k', 'linestyle', 'none', 'linewidth', 1);
errorbar(x(2,:), decodingPowersLower, decodingErrorsLower, 'k', 'linestyle', 'none', 'linewidth', 1);
errorbar(x(3,:), waitingPowersLower, waitingErrorsLower, 'k', 'linestyle', 'none', 'linewidth', 1);

% Customize plot
set(gca, 'XTick', 1:length(channelNamesList));
set(gca, 'XTickLabel', channelNamesList);
legend({'Encoding', 'Decoding', 'Waiting'});
ylabel('Average Gamma Power (Lower Band)');
title('Avg Gamma Power for Encoding v. Decoding v. Waiting (Lower Band)');
hold off;

figure;
barVals = [encodingPowersHigher; decodingPowersHigher; waitingPowersHigher]'; % Grouped bar values
bar(barVals, 'grouped'); % Plot bar chart
hold on;
numGroups = size(barVals, 1);
numBars = size(barVals, 2);
groupWidth = min(0.8, numBars / (numBars + 1.5));
x = nan(numBars, numGroups);
for i = 1:numBars
    x(i,:) = (1:numGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*numBars);
end
errorbar(x(1,:), encodingPowersHigher, encodingErrorsHigher, 'k', 'linestyle', 'none', 'linewidth', 1);
errorbar(x(2,:), decodingPowersHigher, decodingErrorsHigher, 'k', 'linestyle', 'none', 'linewidth', 1);
errorbar(x(3,:), waitingPowersHigher, waitingErrorsHigher, 'k', 'linestyle', 'none', 'linewidth', 1);
set(gca, 'XTick', 1:length(channelNamesList));
set(gca, 'XTickLabel', channelNamesList);
legend({'Encoding', 'Decoding', 'Waiting'});
ylabel('Average Gamma Power (Higher Band)');
title('Avg Gamma Power for Encoding v. Decoding v. Waiting (Higher Band)');
hold off;

%% Helper Functions (same as in MemoryThetaMeasurements.m)
function dataCell = getStateData(tSel, dataOneChannelSel, stateStart, stateEnd)
    dataCell = cell(1, length(stateStart)); % Store each segment separately
    for i = 1:length(stateStart)
        idx = (tSel >= stateStart(i)) & (tSel <= stateEnd(i)); 
        segment = dataOneChannelSel(idx);
        if ~isempty(segment)
            dataCell{i} = segment(:); % Store each segment as a separate cell
        end
    end
end

function [avgPower, powerSEM] = computeGammaPower(dataSegments, ~)
    powerVals = []; % Store power values across all data points
    
    % Collect power values from all segments
    for i = 1:length(dataSegments)
        if ~isempty(dataSegments{i})
            powerVals = [powerVals; dataSegments{i}.^2]; % Square amplitude for power
        end
    end
    
    % Compute mean power
    avgPower = mean(powerVals, 'omitnan'); 
    
    % Compute SEM using total number of data points
    N = length(powerVals); % Total number of data points across segments
    if N > 1
        powerSEM = std(powerVals, 'omitnan') / sqrt(N); % Standard Error of the Mean (SEM)
    else
        powerSEM = 0; % If only one data point, SEM is zero
    end
end

function [stateInd, stateStartTime, stateEndTime] = ...
    findExpState(expState, expStates, SrlTimes, tRel)

% find timing of state start and end 
stateOn = strcmp(expStates, expState); % true when on, false when off 
stateChange = [0; diff(stateOn)]; 
stateEnd = stateChange < 0; % true -> false (-1) means end
stateStart = stateChange > 0; % false -> true (+1) means start 
stateStartTime = SrlTimes(stateStart); stateEndTime = SrlTimes(stateEnd);

% match start and end times 
stateStartTime = sort(stateStartTime); stateEndTime = sort(stateEndTime);
timeEnd = nan(size(stateStartTime)); 
for ti = 1:length(stateStartTime)
    ts = stateStartTime(ti); 
    te = stateEndTime >= ts; 
    te = stateEndTime(te);
    if ~isempty(te)
        te = min(te); 
        if ti < length(stateStartTime)
            if te <= stateStartTime(ti+1)
                timeEnd(ti) = te;
            end
        else
            timeEnd(ti) = te;
        end
    end
end
stateStartTime = stateStartTime(~isnan(timeEnd)); 
stateEndTime = timeEnd(~isnan(timeEnd)); 

% get timing of this state 
stateInd = false(size(tRel)); 
for ti = 1:length(stateStartTime)
    ts = stateStartTime(ti); te = stateEndTime(ti);
    stateInd((tRel > ts) & (tRel < te)) = true;
end

end