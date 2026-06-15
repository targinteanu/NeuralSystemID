%% setup 

% band freq bounds (Hz)
freq1 = [3, 5]; % low theta
%freq1 = [5, 6]; % mid theta
freq2 = [30, 80]; % low gamma
%freq2 = [80, 120]; % high gamma
%freq2 = [120, 180]; % very high gamma

% get folder 
fp = uigetdir; 
FF = dir(fp);

% inst outputs 
tblElec = [];
tblRow = [];

for F = FF'
%% obtain segmented, artifact-free data 

% file selection
%{
[fn,fp] = uigetfile('*SegmentData*.mat', 'Choose Artifact-Free Data File');
load(fullfile(fp,fn));
%}
if ~F.isdir
fn = F.name;
[~,fn,fe] = fileparts(fn);
if strcmpi(fe, '.mat') && contains(fn,'_SegmentData')
load(fullfile(F.folder, F.name));
subjID = fn(1:8);
ArtRemoveDone = contains(fn, '_ArtifactRemoveOffline');
if ArtRemoveDone
    fnOrig = split(fn, '_ArtifactRemoveOffline'); fnOrig = fnOrig{1};
    fpOrig = fp;
    while ~exist(fullfile(fpOrig,[fnOrig,fe]), 'file')
        fpOrig = fileparts(fpOrig); % try one folder out
    end
    load(fullfile(fpOrig,[fnOrig,fe]));
end

tblBL = tblsBaseline{1,1}; tblBL.Properties.Description = 'Baseline';
tblsToAnalyze = [{tblBL}; tblsMisc(:,1)];
tblNames = cellfun(@(T) T.Properties.Description, tblsToAnalyze, 'UniformOutput',false);
tblsSel = listdlg("PromptString","Select for Analysis", "SelectionMode","multiple", ...
    "ListString",tblNames);
tblsToAnalyze = tblsToAnalyze(tblsSel);

%%
for Ti = 1:height(tblsToAnalyze)
%% data processing

% clean baseline selection 
PDdata = tblsToAnalyze{Ti};
%RecordingDescription = PDdata.Properties.Description
srate = PDdata.Properties.SampleRate;
if isnan(srate)
    srate = 1/median(seconds(diff(PDdata.Time)));
end

if Ti > 1
    tblName = PDdata.Properties.Description;
else
    tblName = 'Baseline';
end
tblName = [subjID,'_',tblName];
disp(tblName);

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

% notch out power line noise and any harmonics in band
f0 = 60; % power line fundamental 
qFactor = 35;
notchB = 1; notchA = 1;
for h = f0:f0:min((srate/2), max([freq1,freq2]))
    % add a notch at harmonic h
    [notchBh,notchAh] = iirnotch(h/(srate/2), (h/(srate/2))/qFactor);
    notchB = conv(notchB, notchBh); notchA = conv(notchA, notchAh);
end
figure; freqz(notchB,notchA,[],srate); sgtitle('Power Line Notch Filter');

% filter band of interest 
filtwts = fir1(1024, freq1./(srate/2));
figure; freqz(filtwts,1,[],srate); sgtitle('Band Filter');

pause(.001); drawnow; pause(.001);

% signal processing 
%thetaPowerCortex = bandpower(PDdata.Variables, srate, freqbnd);
%thetaPowerCortex(isoutlier(sqrt(thetaPowerCortex))) = nan;
%PDdata = FilterTimetable(@(b,x) filtfilt(b,1,x), filtwts, PDdata);
for ch = 1:width(PDdata)
    disp(['Filtering channel ',num2str(ch),' of ',num2str(width(PDdata))])
    x = PDdata{:,ch};
    if sum(~isnan(x))
        if (length(notchB) > 1) || (length(notchA) > 1) || (notchB(1) ~= 1)
            x = filtfilt(notchB,notchA,x);
        end
        PDdata{:,ch} = filtfilt(filtwts,1,x);
    end
end
PDdata = PDdata(:,1:63);
%PD_Phase_Data = instPhaseFreqTbl(PDdata);
PD_Channel_Names = PDdata.Properties.VariableNames;

ElecTbl = table('RowNames',PD_Channel_Names(1:63)); % cortical only
RowTbl = table();

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
%thetaPowerWindowed = median(thetaPowerWindowed, 'omitnan');

ElecTbl.([tblName,'_MED']) = median(thetaPowerWindowed, 'omitnan')';
ElecTbl.([tblName,'_STD']) = std(thetaPowerWindowed, 'omitnan')';
ElecTbl.([tblName,'_NUM']) = sum(~isnan(thetaPowerWindowed))';
tblElec = [tblElec, ElecTbl];
%{
G = zeros(size(thetaPowerWindowed,1),21,3);
G(:) = thetaPowerWindowed; 
G = [G(:,:,1); G(:,:,2); G(:,:,3)];
R = mean(G,1, 'omitnan'); % median?
%}
G = zeros(21,3,2);
G(:,:,1) = median(thetaPowerWindowed, 'omitnan')';
G(:,:,2) = std(thetaPowerWindowed, 'omitnan')';
RowTbl.([tblName,'_AVG']) = mean(G(:,:,1),2, 'omitnan'); % ** substitute with wavelet transform based image decomposition
RowTbl.([tblName,'_STD']) = rms(G(:,:,2),2, 'omitnan'); % ** substitute?
RowTbl.([tblName,'_NUM']) = sum(~isnan(G(:,:,1)),2);
tblRow = [tblRow, RowTbl];

OLthresh = thetaPowerWindowed > median(thetaPowerWindowed, 'omitnan');
OLthresh = thetaPowerWindowed(OLthresh); 
io = isoutlier(OLthresh, 'mean'); OLthresh = min(OLthresh(io));
%keyboard; % copy R to separate spreadsheet 

end
end
end

end

% saving: by channels
svname = 'AllChans.xlsx';
svname = inputdlg('Save Channels Data As:', 'Save Channels?', 1, {svname});
if ~isempty(svname)
    svname = svname{1}; 
    writetable(tblElec, fullfile(fp,svname));
end

% saving: by rows
svname = 'Rows.xlsx';
svname = inputdlg('Save Rows Data As:', 'Save Rows?', 1, {svname});
if ~isempty(svname)
    svname = svname{1}; 
    writetable(tblRow, fullfile(fp,svname));
end