%% Calculate and display Theta measurements during Memory task 
% Code by Jack Coursen & Toren Arginteanu

thetaPowerResults = struct();

%% Loading the data 
% 
% At the end of this section, you should have the LFP channels in a
% timetable NS2 sampled at 1kHz, and the corresponding struct ns2. You
% should have the Behnke-Fried microelectrode raw rec in the table NS5
% sampled at 30kHz and corresponding struct ns5. 
% 
% You should also have an event table and structure NEV/nev. This has the
% serial events and the results of the Blackrock machine's spike sorting. I
% am not sure how good this is; if you look at the waveforms, many of the
% "spikes" just look like noise. So I think we will want to ignore it for
% now and try our own spike sorting. 
% 
% Finally, you should have a log of the serial messages "SerialLog"
% obtained through the OnlineDisplaySavedData file. These labels are more
% useful than the ones from the eventdata, but the timing may be less
% precise. Let's use this for now.
% 

filepath = 'Saved Data 2025-05-22 11.44.08'; 
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
nev = openNEV(NEVFile, 'nosave', 'nomat', 'uV'); NEV = nev2table(nev);

%% Looping through the channels 
% 
% This code extracts only the theta power and produces a bar graph showing
% theta power in each channel during memory encoding and decoding. Can you
% modify the code to show gamma power in addition to theta? How about
% theta-gamma coupling (PAC)? Please use the included function calcPAC and
% note that theta is the "Phase" band while gamma is the "Amplitude" band.
% 
% The following code only uses LFP-based measurements from the NS2. 
% In order to incorporate spike data, the NS5 data needs to be processed
% through spike sorting; alternatively, we could try using the Blackrock
% spike sorting in the NEV data.
% 

SamplingFreq = ns2.MetaTags.SamplingFreq;
t = NS2.Time';
tRel = seconds(t-t(1));

% Select channels of interest: 
channelNames = NS2.Properties.VariableNames;
channelIndices = contains(channelNames, 'LH'); % Left Hippocampus
channelIndices = find(channelIndices);

% Stimulation times (pulse train): 
StimTrainRec = NS5.AINP1; % Analog Input 1
% downsample to NS2 (by factor of 30)
StimTrainRec = movmean(StimTrainRec, round(ns5.MetaTags.SamplingFreq/SamplingFreq));
StimTrainRec = interp1(NS5.Time, StimTrainRec, NS2.Time, "nearest");
StimTrainRec = [false; diff(StimTrainRec) > 2000]';

% setup bandpass filter 
bpf = buildFIRBPF(SamplingFreq, 4, 9); % Theta: 4 to 9 Hz

for idx = 1:length(channelIndices)

channelIndex = channelIndices(idx);
channelName = NS2.Properties.VariableNames{channelIndex}
dataOneChannel = NS2{:,channelIndex}';

dataOneChannelWithArtifact = dataOneChannel; 
t0 = t(1);

%% artifact detection
artExtend = 10; % extend artifact by __ samples 
%if numel(channelIndexStim)
    artIndAll = StimTrainRec; % cerestim trigs
    % no other sources of artifact in the memory protocol
%else
%    warning('Stimulus channel ainp1 was not connected.')
%    % assume there are cerestim trigs, but they are not recorded
%    artIndAll = isoutlier(dataOneChannel, 'mean');
%end 
%artIndAll(StimInd) = true;
artIndAll = movsum(artIndAll, artExtend) > 0;
artIndAll_PulseTrain = artIndAll;
artIndAll = find(artIndAll);
[~,baselineStartInd] = max(diff(artIndAll));
baselineEndInd = artIndAll(baselineStartInd+1); baselineStartInd = artIndAll(baselineStartInd); 

%% set baseline to fit model 
if ~isempty(artIndAll)
baselineWinLen = 1000; ARlen = 10; % samples 
dataBaseline = dataOneChannel(baselineStartInd:baselineEndInd); 
dataBaseline = filtfilt(bpf,1,dataBaseline);
baselineWin = (baselineEndInd-baselineStartInd) + [-1,1]*baselineWinLen; 
baselineWin = baselineWin/2; baselineWin = round(baselineWin); 
baselineWin(1) = max(1,baselineWin(1)); baselineWin(2) = min(length(dataBaseline),baselineWin(2));
dataBaseline = dataBaseline(baselineWin(1):baselineWin(2));
ARmdl = ar(iddata(dataBaseline', [], 1/SamplingFreq), ARlen, 'yw');

%% remove artifact 
dataOneChannel = dataOneChannelWithArtifact;
dataOneChannel = dataOneChannel - mean(dataOneChannel); % correct DC offset

for ind = artIndAll
    ind0 = ind - ARlen;
    if ind0 > 0
        dataOneChannel(ind) = myFastForecastAR(ARmdl, dataOneChannel(ind0:(ind-1))', 1);
    end
end
end

% filter 
dataOneChannel = filtfilt(bpf,1,dataOneChannel);

%% select time of interest (manually)
% TO DO: make this automatic, pulled from notes.txt ?
% selind = true(size(t)); % no selection
% select only times before stimulation started: 
iFirstStim = find(StimTrainRec); iFirstStim = iFirstStim(1);
selind = t <= t(iFirstStim) - seconds(60);
tSel = t(selind); tRelSel = tRel(selind);
dataOneChannelSel = dataOneChannel(selind);
selind = find(selind); 

%% determine encode/decode phases of experiment 
expStates = {SerialLog.ParadigmPhase}';
SrlTimes = [SerialLog.TimeStamp]';
SrlTimesSel = (SrlTimes <= max(tRelSel)) & (SrlTimes >= min(tRelSel));
expStates = expStates(SrlTimesSel);
SrlTimesSel = SrlTimes(SrlTimesSel);
[indEncode, encodeStart, encodeEnd] = ...
    findExpState('ENCODE', expStates, SrlTimesSel, tRelSel);
[indDecode, decodeStart, decodeEnd] = ...
    findExpState('DECODE', expStates, SrlTimesSel, tRelSel);
indNeither = (~indEncode)&(~indDecode);

% convert times to absolute 
varnames = {'encodeStart', 'encodeEnd', 'decodeStart', 'decodeEnd'};
for V = varnames
    v = V{:};
    eval([v,' = seconds(',v,') + t0;'])
end

%% Plot time series 
% Individual channel plots have been commented out for now. Alternatively,
% consider making one or two large figures with subplots/tiles for each
% individual channel. 

%{

% select the data to plot 
dataMinMax = [min(dataOneChannelSel), max(dataOneChannelSel)];
dataMinMax = dataMinMax + [-1,1]*.01*diff(dataMinMax);
dataMin = dataMinMax(1); dataMax = dataMinMax(2); 

% plot data  
figure; 
plot(tSel,dataOneChannelSel);
grid on; hold on; lgd = ["data"];

% shade plot regions indicating encode and decode state of paradigm 
if numel(encodeStart)
patch([encodeStart, encodeEnd, encodeEnd, encodeStart]', ...
    repmat([dataMax; dataMax; dataMin; dataMin],1,length(encodeStart)), ...
    'c', 'FaceAlpha', .2); 
lgd = [lgd, "Encode"];
end
if numel(decodeStart)
patch([decodeStart, decodeEnd, decodeEnd, decodeStart]', ...
    repmat([dataMax; dataMax; dataMin; dataMin],1,length(decodeStart)), ...
    'g', 'FaceAlpha', .2); 
lgd = [lgd, "Decode"];
end

% label the plot 
legend(lgd)

%}

%% analyze measurements during encode/decode periods 

encodeData = getStateData(tSel, dataOneChannelSel, encodeStart, encodeEnd);
decodeData = getStateData(tSel, dataOneChannelSel, decodeStart, decodeEnd);

% Compute Theta Power for Encoding and Decoding
[avgThetaPowerEncoding, stdThetaPowerEncoding] = computeThetaPower(encodeData, SamplingFreq);
[avgThetaPowerDecoding, stdThetaPowerDecoding] = computeThetaPower(decodeData, SamplingFreq);

disp(['Avg Theta Power Encoding: ', num2str(avgThetaPowerEncoding), ' ± ', num2str(stdThetaPowerEncoding)]);
disp(['Avg Theta Power Decoding: ', num2str(avgThetaPowerDecoding), ' ± ', num2str(stdThetaPowerDecoding)]);



% Concatenate all encoding segments back-to-back
concatenatedEncodeData = [];
for i = 1:length(encodeData)
    concatenatedEncodeData = [concatenatedEncodeData; encodeData{i}; NaN]; % Add NaN to separate segments
end

% Concatenate all decoding segments back-to-back
concatenatedDecodeData = [];
for i = 1:length(decodeData)
    concatenatedDecodeData = [concatenatedDecodeData; decodeData{i}; NaN]; % Add NaN to separate segments
end

% Find the global min and max for the y-axis range
globalMin = min([concatenatedEncodeData; concatenatedDecodeData], [], 'omitnan');
globalMax = max([concatenatedEncodeData; concatenatedDecodeData], [], 'omitnan');

%{

% Create a figure for encoding and decoding data
figure;
tiledlayout(3,1);

% Plot Encoding Data
nexttile;
plot(concatenatedEncodeData, 'b'); % Blue for encoding
title('Encoding Signal');
ylabel('Amplitude');
xlabel('Time Points');
ylim([globalMin globalMax]); % Set same y-axis range

% Plot Decoding Data
nexttile;
plot(concatenatedDecodeData, 'g'); % Green for decoding
title('Decoding Signal');
ylabel('Amplitude');
xlabel('Time Points');
ylim([globalMin globalMax]); % Set same y-axis range

% Plot Decoding and Encoding Data
nexttile;
plot(concatenatedDecodeData, 'g');
title('Encoding and Decoding Signal (Overlay)');
ylabel('Amplitude');
xlabel('Time Points');
hold on;
plot(concatenatedEncodeData, 'b');

ylim([globalMin globalMax]); % Set same y-axis range

%}

thetaPowerResults(idx).channelName = channelNames{channelIndex}; % Store channel name
thetaPowerResults(idx).encodingPower = avgThetaPowerEncoding;
thetaPowerResults(idx).decodingPower = avgThetaPowerDecoding;
thetaPowerResults(idx).encodingError = stdThetaPowerEncoding;
thetaPowerResults(idx).decodingError = stdThetaPowerDecoding;

end

channelNamesList = {thetaPowerResults.channelName};
encodingPowers = [thetaPowerResults.encodingPower]; 
decodingPowers = [thetaPowerResults.decodingPower]; 
encodingErrors = [thetaPowerResults.encodingError]; 
decodingErrors = [thetaPowerResults.decodingError]; 
%%
figure;
barVals = [encodingPowers; decodingPowers]'; % Grouped bar values
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

% Customize plot
set(gca, 'XTick', 1:length(channelNamesList));
set(gca, 'XTickLabel', channelNamesList);
legend({'Encoding', 'Decoding'});
ylabel('Average Theta Power');
title('Avg Theta Power for Encoding v. Decoding');
hold off;















%% helper functions 


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


function [avgPower, powerSEM] = computeThetaPower(dataSegments, ~)
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
