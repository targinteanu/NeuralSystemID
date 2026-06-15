%% setup 

% choose band freq bounds (Hz)
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
%%
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
PDdata = PDdata(:,1:63); % cortical only 
X = PDdata.Variables; Y = nan([size(X), 2]);
X = X - mean(X,2,'omitnan'); % common avg reref 

% notch out power line noise and any harmonics in band
f0 = 60; % power line fundamental 
qFactor = 35;
notchB = 1; notchA = 1;
for h = f0:f0:min((srate/2), max([freq1,freq2]))
    % add a notch at harmonic h
    [notchBh,notchAh] = iirnotch(h/(srate/2), (h/(srate/2))/qFactor);
    notchB = conv(notchB, notchBh); notchA = conv(notchA, notchAh);
end
%figure; freqz(notchB,notchA,[],srate); sgtitle('Power Line Notch Filter');

% filter bands of interest 
filt1 = fir1(1023, freq1./(srate/2));
filt2 = fir1(1023, freq2./(srate/2));
%figure; freqz(filtwts,1,[],srate); sgtitle('Band Filter');

%pause(.001); drawnow; pause(.001);

% signal processing 
%thetaPowerCortex = bandpower(PDdata.Variables, srate, freqbnd);
%thetaPowerCortex(isoutlier(sqrt(thetaPowerCortex))) = nan;
%PDdata = FilterTimetable(@(b,x) filtfilt(b,1,x), filtwts, PDdata);
for ch = 1:width(PDdata)
    disp(['Filtering channel ',num2str(ch),' of ',num2str(width(PDdata))])
    x = X(:,ch);
    if sum(~isnan(x))
        if (length(notchB) > 1) || (length(notchA) > 1) || (notchB(1) ~= 1)
            x = filtfilt(notchB,notchA,x);
        end
        y1 = filtfilt(filt1,1,x);
        y2 = filtfilt(filt2,1,x);
        Y(:,ch,1) = y1; Y(:,ch,2) = y2; X(:,ch) = x;
    end
end
%PDdata = PDdata(:,1:63);
%PD_Phase_Data = instPhaseFreqTbl(PDdata);
PD_Channel_Names = PDdata.Properties.VariableNames;

ElecTbl = table('RowNames',PD_Channel_Names(1:63)); % cortical only
RowTbl = table();

%% windowed power 

% Calculate windowed power
windowSize = ceil(1 * srate); % samples 
powerWindowed = Y.^2;
for ch = 1:size(Y,2) % channel
    disp(['Calculating power: channel ',num2str(ch),' of ',num2str(width(PDdata))])
    for f = 1:size(Y,3) % freq band (i.e. theta, gamma)
        x = powerWindowed(:,ch,f);
        if sum(~isnan(x))
            powerWindowed(:,ch,f) = movmean(envelope(x.^2), windowSize);
        end
    end
end
powerWindowed = powerWindowed(round(windowSize/2):windowSize:end, :, :);
powerWindowed = 10*log10((powerWindowed)); % decibel (dB) scale 
powerWindowed(isoutlier(powerWindowed)) = nan;
%powerWindowed = median(powerWindowed, 'omitnan');
thetaPowerWindowed = powerWindowed(:,:,1); 
gammaPowerWindowed = powerWindowed(:,:,2);

% store THETA
ElecTbl.([tblName,'_The_MED']) = median(thetaPowerWindowed, 'omitnan')';
ElecTbl.([tblName,'_The_STD']) = std(thetaPowerWindowed, 'omitnan')';
ElecTbl.([tblName,'_The_NUM']) = sum(~isnan(thetaPowerWindowed))';
G = makegrid(thetaPowerWindowed);
RowTbl.([tblName,'_The_AVG']) = mean(G(:,:,1),2, 'omitnan'); % ** substitute with wavelet transform based image decomposition
RowTbl.([tblName,'_The_STD']) = rms(G(:,:,2),2, 'omitnan'); % ** substitute?
RowTbl.([tblName,'_The_NUM']) = sum(~isnan(G(:,:,1)),2);

% store GAMMA 
ElecTbl.([tblName,'_Gam_MED']) = median(gammaPowerWindowed, 'omitnan')';
ElecTbl.([tblName,'_Gam_STD']) = std(gammaPowerWindowed, 'omitnan')';
ElecTbl.([tblName,'_Gam_NUM']) = sum(~isnan(gammaPowerWindowed))';
G = makegrid(gammaPowerWindowed);
RowTbl.([tblName,'_Gam_AVG']) = mean(G(:,:,1),2, 'omitnan'); % ** substitute with wavelet transform based image decomposition
RowTbl.([tblName,'_Gam_STD']) = rms(G(:,:,2),2, 'omitnan'); % ** substitute?
RowTbl.([tblName,'_Gam_NUM']) = sum(~isnan(G(:,:,1)),2);

%{
OLthresh = thetaPowerWindowed > median(thetaPowerWindowed, 'omitnan');
OLthresh = thetaPowerWindowed(OLthresh); 
io = isoutlier(OLthresh, 'mean'); OLthresh = min(OLthresh(io));
%keyboard; % copy R to separate spreadsheet 
%}

%% windowed PAC 

% hilbert 
Z = Y; 
Z(:,:,1) = angle(hilbert(Z(:,:,1)));
Z(:,:,2) =   abs(hilbert(Z(:,:,2)));

% windowed calc 
winStart = 1:windowSize:size(Z,1); winEnd = winStart+windowSize;
while winEnd(end) > size(Z,1) % trim incomplete windows 
    winStart = winStart(1:(end-1)); winEnd = winEnd(1:(end-1));
end
PACwindowed = nan(length(winStart),size(Z,2));
for ch = 1:size(Z,2) % channel
    disp(['Calculating PAC: channel ',num2str(ch),' of ',num2str(width(PDdata))])
    for w = 1:length(winStart)
        phi = Z(winStart(w):winEnd(w),ch,1); 
        amp = Z(winStart(w):winEnd(w),ch,2); 
        PACwindowed(w,ch) = calcPAChelper(phi,amp,18); % 18 bins
    end
end

% exclude some windows to avoid edge effects 
trimwin = 0.1*size(Z,1); % at most central 80% 
trimwin = ceil(trimwin/windowSize);
PACwindowed = PACwindowed((trimwin+1):(end-trimwin),:);

% store PAC 
ElecTbl.([tblName,'_PAC_MED']) = median(PACwindowed, 'omitnan')';
ElecTbl.([tblName,'_PAC_STD']) = std(PACwindowed, 'omitnan')';
ElecTbl.([tblName,'_PAC_NUM']) = sum(~isnan(PACwindowed))';
G = makegrid(PACwindowed);
RowTbl.([tblName,'_PAC_AVG']) = mean(G(:,:,1),2, 'omitnan'); % ** substitute with wavelet transform based image decomposition
RowTbl.([tblName,'_PAC_STD']) = rms(G(:,:,2),2, 'omitnan'); % ** substitute?
RowTbl.([tblName,'_PAC_NUM']) = sum(~isnan(G(:,:,1)),2);

%%
tblElec = [tblElec, ElecTbl];
tblRow = [tblRow, RowTbl];

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

%% helper(s) 

function G = makegrid(allchandata)
%{
G = zeros(size(allchandata,1),21,3);
G(:) = allchandata; 
G = [G(:,:,1); G(:,:,2); G(:,:,3)];
R = mean(G,1, 'omitnan'); % median?
%}
G = zeros(21,3,2);
g = zeros(21,3);
g(:) = median(allchandata, 'omitnan')';
G(:,:,1) = g;
g(:) = std(allchandata, 'omitnan')';
G(:,:,2) = g;
end