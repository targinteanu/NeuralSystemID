% sensitivity analysis comparing dynamic AR parameters 

% file selection
thisfilename = mfilename("fullpath");
[fn,fp] = uigetfile('*SegmentData*.mat', 'Choose Segmented Data File');
SegmentedDataFullfile = fullfile(fp,fn);
load(SegmentedDataFullfile);
[~,fn,fe] = fileparts(fn);
subjname = upper(fn(1:8));
isEMU = contains(subjname, 'PY');

% data selection 
tblBaseline = tblsBaseline{1};
chName = tblBaseline.Properties.VariableNames; 
chDesc = tblBaseline.Properties.VariableDescriptions;
chNameOnly = regexprep(upper(chName), '-REF', '');
chNameOnly = regexprep(chNameOnly, '\d+', ''); % no number
if isEMU
% only include macro channels both named (i.e. targeted) and identified for
% hippocampus with no uncertainty
chsel = ...
    (strcmp(chNameOnly,'LH')  | strcmp(chNameOnly,'RH')  | ... R/L hippo
     strcmp(chNameOnly,'LAH') | strcmp(chNameOnly,'RAH') | ... R/L ant hippo
     strcmp(chNameOnly,'LPH') | strcmp(chNameOnly,'RPH')) & ... R/L post hippo
    contains(lower(chDesc),'hip') & ~contains(chDesc,'?'); ...
%    & ~contains(chName,'BF'); 
else
    chsel = contains(chDesc, 'PrG'); % Precentral Gyrus
end
dtaBL = tblBaseline(:,chsel);
chselName = dtaBL.Properties.VariableNames;

% sample rate 
Fs = dtaBL.Properties.SampleRate;
if isnan(Fs)
    dt = seconds(diff(dtaBL.Time));
    if (min(dt) < 0) 
        error('Time series must be ascending.')
    end
    dtmean = mean(dt);
    if max(abs(dt-dtmean)) > .005
        error('Sample rate is not consistent.')
    end
    Fs = 1/dtmean;
end

%% subselect channels 
% try most and least average channels 
X = dtaBL.Variables; 
X = X - mean(X,1); X = X - mean(X,2); 
X = sum(X.^2, 1);
[~,chsel1] = min(X); [~,chsel2] = max(X);
dtaBL = dtaBL(:, [chsel1, chsel2]);
chselName = dtaBL.Properties.VariableNames;

%% parameters 
packetSize = 20; % samples
phTargets = [0, pi]; % radians
ARwin = 1000; predWin = 500; % samples
ARord = 50; % # coefficients 
if isEMU
    freqrng = [4, 9]; % Hz (Theta band)
else
    freqrng = [13, 30]; % Hz (beta band)
end

learnrates = .05*2.^(-3:3); 

%% loop / parameter sweep 
errResults = cell(length(phTargets), width(dtaBL), length(learnrates));
numResults = cell(size(errResults));
trgResults = cell(size(errResults));
% dim 1: target phase
% dim 2: channel 
% dim 3: learn rate

for c = 1:1:width(dtaBL)

chtoplot = dtaBL.Properties.VariableNames{c}

for p = 1:length(phTargets)

phTarget = phTargets(p)

parfor l = 1:length(learnrates)

learnrate = learnrates(l)

prog = ((c-1)*length(phTargets) + p-1 )*length(learnrates) + l-1;
prog = prog/(length(phTargets)*width(dtaBL)*length(learnrates));
disp([' ========== PROGRESS: ',num2str(round(100*prog)),'% ========== ']);

% run phase detection 
[phAll, phEst, frAll, frEst, trgTimeDyn, phStimDyn, ~, ~, durD, ...
    nCycleDyn, nMissDyn, nXtraDyn] = ...
    offline_PhaseDetect(dtaBL.(chtoplot)', Fs, [], dtaBL.Time', chtoplot, ...
    phTarget, freqrng, ARwin, ARord, predWin, -1, packetSize, ...
    learnrate, false, false, false, false, false);
phErrDyn = radfix(phEst-phAll); frErrDyn = frEst - frAll;

pause(.001); 
drawnow;
pause(.001);

errResults{p,c,l} = radfix(phStimDyn-phTarget);
trgResults{p,c,l} = trgTimeDyn;
numResults{p,c,l} = [nMissDyn, nXtraDyn, nCycleDyn];

end % learn rates 
end % phase targets
end % channels

%% aggregate/display all channel/target results 

errResultsAll = cell(1, size(errResults,3));
numResultsAll = cell(1, size(numResults,3));
for r = 1:size(errResults,1)
    for h = 1:size(errResults,3)
        for c = 1:size(errResults,2)
            errResultsAll{1,h} = [errResultsAll{1,h}, errResults{r,c,h}];
            numResultsAll{1,h} = [numResultsAll{1,h}; numResults{r,c,h}];
        end
    end
end
for h = 1:size(numResults,3)
    numResultsAll{1,h} = sum(numResultsAll{1,h});
    numResultsAll{1,h} = [numResultsAll{1,h}(1:2), ...
        numResultsAll{1,h}(3)-sum(numResultsAll{1,h}(1:2))];
end

errResultsAll_avg = cellfun(@(ph) circ_mean(ph'), errResultsAll);
errResultsAll_std = cellfun(@(ph) circ_std(ph'),  errResultsAll);

fig1 = figure; 
for h = 1:size(errResultsAll,2)
    polarhistogram(errResultsAll{h}, 18); hold on; 
end
title('Phase Error (causal - offline)'); 
lgd = legend(string(learnrates)); lgd.Title.String = 'Learn Rate';

fig2 = figure; 
bar(learnrates, errResultsAll_avg); 
hold on; grid on; 
errorbar(learnrates,errResultsAll_avg, errResultsAll_std,errResultsAll_std, '.');
xlabel('Learn Rate'); ylabel('Phase Error (rad)'); 
legend('Circ Mean', '±1 Circ SD');

%% save 
thisfilever = getFileVersion(thisfilename);
[~,thisfilename] = fileparts(thisfilename);
svname = [fp,fn,'_',thisfilename,'_',thisfilever];

save(svname, 'errResults', 'numResults', 'trgResults', ...
    'phTargets', 'chselName', 'packetSize', 'ARwin', 'ARord', ...
    'learnrates', 'freqrng');
saveas(fig1, [svname,'_polhisto'], 'png');
saveas(fig2, [svname,'_rectbar'], 'png');