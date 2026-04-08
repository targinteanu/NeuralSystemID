%% obtain segmented, artifact-free data 

% file selection
[fn,fp] = uigetfile('*SegmentData*.mat', 'Choose Artifact-Free Data File');
load(fullfile(fp,fn));
[~,fn,fe] = fileparts(fn);
ArtRemoveDone = contains(fn, '_ArtifactRemoveOffline');
if ArtRemoveDone
    fnOrig = split(fn, '_ArtifactRemoveOffline'); fnOrig = fnOrig{1};
    fpOrig = fp;
    while ~exist(fullfile(fpOrig,[fnOrig,fe]), 'file')
        fpOrig = fileparts(fpOrig); % try one folder out
    end
    load(fullfile(fpOrig,[fnOrig,fe]));
end

tblsToAnalyze = [tblsBaseline(1,1); tblsMisc(:,1)];
tblOut = [];

%%
for Ti = 1:height(tblsToAnalyze)
%% data processing

% clean baseline selection 
PDdata = tblsToAnalyze{Ti};
srate = PDdata.Properties.SampleRate;
if isnan(srate)
    srate = 1/median(seconds(diff(PDdata.Time)));
end

if Ti > 1
    tblName = PDdata.Properties.Description;
else
    tblName = 'Baseline';
end

%% channel selection
channames = num2str((1:63)');
channames(channames==' ') = '0';
channames = "LS"+string(channames);
channames = channames';
for ch = channames
    if ~sum(strcmp(PDdata.Properties.VariableNames, ch))
        % fill missing channel with nan
        PDdata.(ch) = nan(height(PDdata),1);
    end
end
PDdata = PDdata(:,channames); % reorder

% common avg ECoG reref 
PDdata{:,:} = PDdata{:,:} - mean(PDdata{:,1:63}, 2, 'omitnan');

% signal processing 
%thetaPowerCortex = bandpower(PDdata.Variables, srate, [4,9]);
%thetaPowerCortex(isoutlier(sqrt(thetaPowerCortex))) = nan;
filtwts = fir1(1024, [4, 9]./(srate/2));
%PDdata = FilterTimetable(@(b,x) filtfilt(b,1,x), filtwts, PDdata);
for ch = 1:width(PDdata)
    disp(['Filtering channel ',num2str(ch),' of ',num2str(width(PDdata))])
    x = PDdata{:,ch};
    if sum(~isnan(x))
        PDdata{:,ch} = filtfilt(filtwts,1,x);
    end
end
PDdata = PDdata(:,1:63);
%PD_Phase_Data = instPhaseFreqTbl(PDdata);
PD_Channel_Names = PDdata.Properties.VariableNames;

ElecTbl = table('RowNames',PD_Channel_Names(1:63)); % cortical only

%% windowed theta power 
% Calculate windowed theta power
windowSize = 1000; % samples 
thetaPowerWindowed = PDdata.Variables.^2;
for ch = 1:width(PDdata)
    disp(['Calculating power: channel ',num2str(ch),' of ',num2str(width(PDdata))])
    x = thetaPowerWindowed(:,ch);
    if sum(~isnan(x))
        thetaPowerWindowed(:,ch) = movmean(envelope(x.^2), windowSize);
    end
end
%thetaPowerWindowed = movmean(envelope(PDdata.Variables).^2, windowSize);
thetaPowerWindowed = thetaPowerWindowed(round(windowSize/2):windowSize:end, :);
thetaPowerWindowed = 10*log10((thetaPowerWindowed)); % decibel (dB) scale 
thetaPowerWindowed(isoutlier(thetaPowerWindowed)) = nan;
thetaPowerWindowed = median(thetaPowerWindowed, 'omitnan');

ElecTbl.(tblName) = thetaPowerWindowed';
tblOut = [tblOut, ElecTbl];
G = zeros(21,3);
G(:) = thetaPowerWindowed; 
R = mean(G,2, 'omitnan');

OLthresh = thetaPowerWindowed > median(thetaPowerWindowed, 'omitnan');
OLthresh = thetaPowerWindowed(OLthresh); 
io = isoutlier(OLthresh, 'mean'); OLthresh = min(OLthresh(io));

end

writetable(tblOut, fullfile(fp,'ThetaPowerAllChans.xlsx'));