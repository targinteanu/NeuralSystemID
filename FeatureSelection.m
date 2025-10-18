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
bandbounds = [0.5,4,9,13,30,70,150];
bandnames = ["\delta", "\theta", "\alpha", "\beta", "lo\gamma", "hi\gamma"];

SampleRate = SampleRates(1);
[pBL, fBL] = periodogram(dtaBL.Variables, [], [], SampleRate, 'power');
pBL = 10*log10(pBL);
figure; semilogx(fBL, pBL); grid on;
xlim([0 200]);
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
title('baseline PSD power estimate');
xticks(bandbounds)

%% organize into tables 

featnames = ["Raw", "Filtered", "Hilbert", "MagPhase"];

alltbls = cell(4,length(featnames));
% stim type {baseline, other non-stim, cort stim, depth stim} x feature
% all wrapped in cells so there can be multiple 
alltbls{1,1} = {dtaBL};

% collect by stim type 
tblsToOrganize = [...
    tblsMisc(:,1); ... selected channels only!!!
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

% TO DO: if any of the above tables has missing or nan for a time period
% greater than tol, break it into cell array of tables 

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
    disp(['Notch Filtering: ',num2str(Ti),' of ',num2str(height(alltbls))])
    for Tj = 1:height(alltbls{Ti,1})
        T = alltbls{Ti,1}{Tj};
        for c = 1:width(T)
            T{:,c} = filtfilt(notchB,notchA,T{:,c});
        end
        alltbls{Ti,1}{Tj} = T;
    end
end

%%{
% reref to common average 
for Ti = 1:height(alltbls)
    disp(['Rereferencing: ',num2str(Ti),' of ',num2str(height(alltbls))])
    for Tj = 1:height(alltbls{Ti,1})
        T = alltbls{Ti,1}{Tj};
        T.Variables = T.Variables - mean(T.Variables, 2, 'omitnan');
        alltbls{Ti,1}{Tj} = T;
    end
end
%}

% visualize results 
sigBL = alltbls{1,1}{1};
figure; periodogram(sigBL.Variables, [], [], SampleRate, 'power');

%% calc features 

% filter 
for Ti = 1:height(alltbls)
    disp(['Band Filtering: ',num2str(Ti),' of ',num2str(height(alltbls))])
    TTraw = alltbls{Ti,1};
    TTfilt = cell(size(TTraw));
    for Tj = 1:height(TTraw)
        Traw = TTraw{Tj};
        Tfilt = [];
        for b = 1:length(bandnames)
            bpf = buildFIRBPF(SampleRate, bandbounds(b), bandbounds(b+1));
            varnames = string(Traw.Properties.VariableNames)+" "+bandnames(b);
            Xfilt = filtfilt(bpf,1,Traw.Variables);
            Tfilt = [Tfilt, array2timetable(Xfilt,...
                "RowTimes",Traw.Time, "VariableNames",varnames)];
        end
        TTfilt{Tj} = Tfilt;
    end
    alltbls{Ti,2} = TTfilt;
end

% hilbert, mag/phase
for Ti = 1:height(alltbls)
    disp(['Hilbert: ',num2str(Ti),' of ',num2str(height(alltbls))])
    TTfilt = alltbls{Ti,2};
    TThilb = cell(size(TTfilt)); TTmaph = cell(size(TTfilt));
    for Tj = 1:height(TTfilt)
        Tfilt = TTfilt{Tj};
        varnames = string(Tfilt.Properties.VariableNames);
        Xhilb = hilbert(Tfilt.Variables);
        TThilb{Tj} = [...
            array2timetable([real(Xhilb)], "RowTimes",Tfilt.Time, ...
                "VariableNames",varnames+" real"), ...
            array2timetable([imag(Xhilb)], "RowTimes",Tfilt.Time, ...
                "VariableNames",varnames+" imag")];
        TTmaph{Tj} = [...
            array2timetable([abs(Xhilb)], "RowTimes",Tfilt.Time, ...
                "VariableNames",varnames+" mag"), ...
            array2timetable([angle(Xhilb)], "RowTimes",Tfilt.Time, ...
                "VariableNames",varnames+" phase")];
    end
    alltbls{Ti,3} = TThilb; alltbls{Ti,4} = TTmaph;
end

%% evaluate difference between stim and baseline 

figure('Units','normalized', 'Position',[.05,.05,.9,.9]); 
pause(.001); drawnow; pause(.001);
for feat = 1:size(alltbls,2)
    tblsBL = alltbls{1,feat};
    T = [];
    for Tj = 1:height(tblsBL)
        T = [T; tblsBL{Tj}];
    end
    XBL = T.Variables;
    varnames = T.Properties.VariableNames;
    for stim = 2:size(alltbls,1)
        subplot(size(alltbls,1)-1, size(alltbls,2), size(alltbls,2)*(stim-2)+feat);
        tbls = alltbls{stim,feat};
        T = [];
        for Tj = 1:height(tbls)
            T = [T; tbls{Tj}];
        end
        X = T.Variables;
        p = nan(1,width(X));
        for c = 1:width(X)
            [~,p(c)] = kstest2(XBL(:,c), X(:,c));
        end
        stem((p)); grid on; 
        if length(p) > length(chanselidx)
            xticks(1:length(chanselidx):length(p)); 
            xticklabels(varnames(1:length(chanselidx):length(p)));
        else
            xticks(1:length(p)); xticklabels(varnames);
        end
        if stim == 2
            title(featnames(feat));
        end
        if feat == 1
            ylabel({char(T.Properties.Description), 'p vs baseline'});
        end
        if stim == size(alltbls,1)
            xlabel('channel/measurement');
        end
        pause(.001); drawnow; pause(.001);
    end
end

%% organize into regularly spaced arrays for each hzn
% TO DO: disp output progress 

hznsN = ceil(hzns*SampleRate); % samples 

allarr = cell([size(alltbls), length(hzns)]);
for feat = 1:size(allarr,2)
    for stim = 1:size(allarr,1)
        tbls = alltbls{stim,feat};
        for hzn = 1:size(allarr,3)
            disp(['Dataset: ',num2str((feat-1)*size(allarr,2)+stim),' of ',num2str(numel(alltbls)),...
                '; horizon ',num2str(hzn),' of ',num2str(size(allarr,3))])
            inp = []; oup = [];
            for Tj = 1:height(tbls)
                % TO DO: this assumes all tbls columns are same order!
                X = tbls{Tj}.Variables; % assume reg spaced time 
                %for startind = 0:(hznsN(hzn)-1)
                startind = 0;
                    inpoup = X((1:hznsN(hzn):end)+startind, :);
                    inp = [inp; inpoup(1:(end-1),:)];
                    oup = [oup; inpoup(2:end,:)];
                %end
            end
            allarr{stim,feat,hzn} = cat(3,inp,oup);
        end
    end
end

%% evaluate "autonomous" linear relationship of "states" 

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