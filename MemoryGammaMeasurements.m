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

% Setup gamma bandpass filter 
gammaBPF = buildFIRBPF(SamplingFreq, 30, 100); % Gamma: 30 to 100 Hz

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
dataBaseline = filtfilt(gammaBPF,1,dataBaseline);
baselineWin = (baselineEndInd-baselineStartInd) + [-1,1]*baselineWinLen; 
baselineWin = baselineWin/2; baselineWin = round(baselineWin); 
baselineWin(1) = max(1,baselineWin(1)); baselineWin(2) = min(length(dataBaseline),baselineWin(2));
dataBaseline = dataBaseline(baselineWin(1):baselineWin(2));
ARmdl = ar(iddata(dataBaseline', [], 1/SamplingFreq), ARlen, 'yw');

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

% Filter for gamma band
dataOneChannel = filtfilt(gammaBPF,1,dataOneChannel);

%% Select time of interest (manually)
iFirstStim = find(StimTrainRec); iFirstStim = iFirstStim(1);
selind = t <= t(iFirstStim) - seconds(60);
tSel = t(selind); tRelSel = tRel(selind);
dataOneChannelSel = dataOneChannel(selind);
selind = find(selind);

%% Determine encode/decode phases of experiment 
expStates = {SerialLog.ParadigmPhase}';
SrlTimes = [SerialLog.TimeStamp]';
SrlTimesSel = (SrlTimes <= max(tRelSel)) & (SrlTimes >= min(tRelSel));
expStates = expStates(SrlTimesSel);
SrlTimesSel = SrlTimes(SrlTimesSel);
[indEncode, encodeStart, encodeEnd] = ...
    findExpState('ENCODE', expStates, SrlTimesSel, tRelSel);
[indDecode, decodeStart, decodeEnd] = ...
    findExpState('DECODE', expStates, SrlTimesSel, tRelSel);
[indWait, waitStart, waitEnd] = ...
    findExpState('WAIT', expStates, SrlTimesSel, tRelSel);

% convert times to absolute 
varnames = {'encodeStart', 'encodeEnd', 'decodeStart', 'decodeEnd', 'waitStart', 'waitEnd'};
for V = varnames
    v = V{:};
    eval([v,' = seconds(',v,') + t0;'])
end

%% Analyze measurements during encode/decode periods 
encodeData = getStateData(tSel, dataOneChannelSel, encodeStart, encodeEnd);
decodeData = getStateData(tSel, dataOneChannelSel, decodeStart, decodeEnd);
waitData = getStateData(tSel, dataOneChannelSel, waitStart, waitEnd);

% Compute Gamma Power for Encoding and Decoding
[avgGammaPowerEncoding, stdGammaPowerEncoding] = computeGammaPower(encodeData, SamplingFreq); % can do calc pac inside computegammapower
[avgGammaPowerDecoding, stdGammaPowerDecoding] = computeGammaPower(decodeData, SamplingFreq);
[avgGammaPowerWaiting, stdGammaPowerWaiting] = computeGammaPower(waitData, SamplingFreq);

disp(['Avg Gamma Power Encoding: ', num2str(avgGammaPowerEncoding), ' ± ', num2str(stdGammaPowerEncoding)]);
disp(['Avg Gamma Power Decoding: ', num2str(avgGammaPowerDecoding), ' ± ', num2str(stdGammaPowerDecoding)]);
disp(['Avg Gamma Power Waiting: ', num2str(avgGammaPowerWaiting), ' ± ', num2str(stdGammaPowerWaiting)]);

gammaPowerResults(idx).channelName = channelNames{channelIndex}; % Store channel name
gammaPowerResults(idx).encodingPower = avgGammaPowerEncoding;
gammaPowerResults(idx).decodingPower = avgGammaPowerDecoding;
gammaPowerResults(idx).encodingError = stdGammaPowerEncoding;
gammaPowerResults(idx).decodingError = stdGammaPowerDecoding;
gammaPowerResults(idx).filteredSignal = dataOneChannel; % Store filtered signal for PAC
gammaPowerResults(idx).encodingSignal = encodeData;
gammaPowerResults(idx).decodingSignal = decodeData;
gammaPowerResults(idx).waitingPower = avgGammaPowerWaiting;
gammaPowerResults(idx).waitingError = stdGammaPowerWaiting;



end

channelNamesList = {gammaPowerResults.channelName};
encodingPowers = [gammaPowerResults.encodingPower]; 
decodingPowers = [gammaPowerResults.decodingPower]; 
encodingErrors = [gammaPowerResults.encodingError]; 
decodingErrors = [gammaPowerResults.decodingError];
waitingPowers = [gammaPowerResults.waitingPower];
waitingErrors = [gammaPowerResults.waitingError];

% Filter out LH03 from the results
excludeChannel = 'LH03';
includeIndices = ~strcmp(channelNamesList, excludeChannel);

% Apply the filter
channelNamesList = channelNamesList(includeIndices);
encodingPowers = encodingPowers(includeIndices);
decodingPowers = decodingPowers(includeIndices);
encodingErrors = encodingErrors(includeIndices);
decodingErrors = decodingErrors(includeIndices);
waitingPowers = waitingPowers(includeIndices);
waitingErrors = waitingErrors(includeIndices);
%%
figure;
barVals = [encodingPowers; decodingPowers; waitingPowers]'; % Grouped bar values
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
errorbar(x(1,:), encodingPowers, encodingErrors, 'k', 'linestyle', 'none', 'linewidth', 1);
errorbar(x(2,:), decodingPowers, decodingErrors, 'k', 'linestyle', 'none', 'linewidth', 1);
errorbar(x(3,:), waitingPowers, waitingErrors, 'k', 'linestyle', 'none', 'linewidth', 1);

% Customize plot
set(gca, 'XTick', 1:length(channelNamesList));
set(gca, 'XTickLabel', channelNamesList);
legend({'Encoding', 'Decoding', 'Waiting'});
ylabel('Average Gamma Power');
title('Avg Gamma Power for Encoding v. Decoding v. Waiting');
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