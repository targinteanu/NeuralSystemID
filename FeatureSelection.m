%% Feature Selection Tool
% Explore preprocessing steps to see which features depend on (A) previous
% values of the same features and (B) stimulation 

hzns = [.05, .1, .5, 1]; % s

%% obtain segmented, artifact-removed data 

% init empty tbl lists 
tblsTrig_ArtifactRemoved = {};
tblsTrig = {};
tblsStimNoTrig_ArtifactRemoved = {};
tblsStimNoTrig = {};

% file selection
thisfilename = mfilename("fullpath");
[fn,fp] = uigetfile('*SegmentData*.mat', 'Choose Artifact-Free Data File');
load(fullfile(fp,fn));
[~,fn,fe] = fileparts(fn);
ArtRemoveDone = contains(fn, '_ArtifactRemoveOffline');
if ArtRemoveDone
    fnOrig = split(fn, '_ArtifactRemoveOffline'); fnOrig = fnOrig{1};
    load(fullfile(fp,[fnOrig,fe]));
end

% select channels of interest 
dtaBL = tblsBaseline{1,1};
chanlist = dtaBL.Properties.VariableNames; chanlist = string(chanlist);
chandesc = dtaBL.Properties.VariableDescriptions; chandesc = string(chandesc);
if ArtRemoveDone
    chanselidx = arrayfun(@(chname) find(strcmpi(chanlist, chname)), chanlistsel);
else
    chanselidx = 1;
end
[chanselidx,chanselmade] = listdlg("ListString",chanlist+": "+chandesc,...
    "SelectionMode","multiple", "PromptString","Select channels to INCLUDE", ...
    "ListSize",[500,300], "InitialValue",chanselidx);
chanlistsel = chanlist(chanselidx);
if chanselmade
    dtaBL = dtaBL(:,chanselidx);
end

%% define freq bands
SampleRate = SampleRates(1);
[pBL, fBL] = periodogram(dtaBL.Variables, [], [], SampleRate, 'power');
pBL = 10*log10(pBL);
figure; semilogx(fBL, pBL); grid on;
xlim([0 200]);
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
title('baseline PSD power estimate');
xticks([0.5,4,9,13,30,70,150])

%% organize into tables 

featnames = ["Raw", "Filtered", "Hilbert", "MagPhase"];

alltbls = cell(4,length(featnames));
% stim type {baseline, other non-stim, cort stim, depth stim} x feature
% all wrapped in cells so there can be multiple 
alltbls{1,1} = {dtaBL};

% collect by stim type 
tblsToOrganize = [...
    tblsMisc(:,1); ...
    selectTbls(tblsTrig, tblsTrig_ArtifactRemoved, chanselmade, chanlistsel); ...
    selectTbls(tblsStimNoTrig, tblsStimNoTrig_ArtifactRemoved, chanselmade, chanlistsel)];
tblsToOrganizeDescs = cellfun(@(T) string(T.Properties.Description), tblsToOrganize);

% other non-stim 
[selidx,selmade] = listdlg("ListString",tblsToOrganizeDescs, ...
    "SelectionMode","multiple", "ListSize",[500,300], ...
    "PromptString","Select NON-stimulation");
if selmade
    alltbls{2,1} = tblsToOrganize(selidx);
end

% cort stim 
[selidx,selmade] = listdlg("ListString",tblsToOrganizeDescs, ...
    "SelectionMode","multiple", "ListSize",[500,300], ...
    "PromptString","Select CORTICAL stimulation");
if selmade
    alltbls{3,1} = tblsToOrganize(selidx);
end

% depth stim 
[selidx,selmade] = listdlg("ListString",tblsToOrganizeDescs, ...
    "SelectionMode","multiple", "ListSize",[500,300], ...
    "PromptString","Select DEPTH stimulation");
if selmade
    alltbls{4,1} = tblsToOrganize(selidx);
end

%% processing steps for all conditions 

% notch out power line noise and harmonics
f0 = 60; % power line fundamental 
qFactor = 35;
notchB = 1; notchA = 1;
for h = f0:f0:(SampleRate/2)
    % add a notch at harmonic h
    [notchBh,notchAh] = iirnotch(h/(SampleRate/2), (h/(SampleRate/2))/qFactor);
    notchB = conv(notchB, notchBh); notchA = conv(notchA, notchAh);
end
for Ti = 1:height(alltbls)
    for Tj = 1:height(alltbls{Ti,1})
        T = alltbls{Ti,1}{Tj};
        for c = 1:width(T)
            T{:,c} = filtfilt(notchB,notchA,T{:,c});
        end
        alltbls{Ti,1}{Tj} = T;
    end
end

% visualize results 
sigBL = alltbls{1,1}{1};
figure; periodogram(sigBL.Variables, [], [], SampleRate, 'power');

%% calc features 

%% organize into regularly spaced arrays for each hzn

%% helpers 

function tblslist = selectTbls(tblsOrig, tblsArtrem, chanselmade, chanlistsel)
% Select between original and art removed table(s). If selected channels
% have not been artifact removed, this will try to pull the original data,
% in which case output table cols might not be in the same order. 
tblsOrig = tblsOrig(:,1);
if ~isempty(tblsArtrem)
    tblslist = tblsArtrem;
    if chanselmade
        for Ti = 1:length(tblslist)
            T = tblslist{Ti};
            for ch = chanlistsel
                if ~sum(strcmp(T.Properties.VariableNames, ch))
                    % current channel ch is desired but not art removed
                    terr = nan(size(tblsOrig));
                    for Tj = 1:length(tblsOrig)
                        To = tblsOrig{Tj};
                        terr(Tj) = sum(...
                            abs(seconds(min(T.Time)-min(To.Time)) + ...
                            abs(seconds(max(T.Time)-max(To.Time)))));
                    end
                    [~,Tji] = min(terr); % find the original table with most time overlap
                    To = tblsOrig{Tji};
                    To_ch = array2timetable(To.(ch), ...
                        "VariableNames",ch, "RowTimes",To.Time);
                    T = synchronize(T,To_ch,'first','nearest');
                end
            end
            tblslist{Ti} = T;
        end
    end
else
    tblslist = tblsOrig;
end
end