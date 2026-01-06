%% obtain segmented, artifact-removed data 

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

%% select channels of interest 
dtaBL = tblsBaseline{1,1};
dtaBL.Properties.Description = 'Baseline';
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
bandbounds = [0.1,2,6,14,40,80,100];
bandnames = ["\delta", "\theta", "\alpha", "\beta", "lo\gamma", "hi\gamma"];
bandcent = .5 * (bandbounds(2:end) + bandbounds(1:(end-1)));

% calc power spectrum 
SampleRate = SampleRates(1);
[pBL, fBL] = periodogram(dtaBL.Variables, [], [], SampleRate, 'power');
pBL = 10*log10(pBL);

% correct for pink noise 
F = [log10(fBL), ones(size(fBL))];
pinkcoef = F(2:end,:)\pBL(2:end,:); 
    % pink noise: P = k*f^a where a < 0
    % 1: slope of log-log, i.e. exponent a of f
    % 2: intercept of log-log, i.e. log10(k)
P = F*pinkcoef; % best-fit pink noise power 

% plot 
pause(.5); fig_spec1 = figure('Units','normalized', 'Position',[.05,.05,.5,.9]); 
subplot(2,1,1); semilogx(fBL, pBL); grid on;
xlim([0 200]);
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
title('baseline PSD power estimate');
xticks(bandbounds)
subplot(2,1,2); semilogx(fBL, pBL-P); grid on;
xlim([0 200]);
xlabel('Frequency (Hz)'); ylabel('Power Diff (dB)');
title('pink noise corrected');
xticks(bandbounds)
pause(.001); drawnow; pause(.001);