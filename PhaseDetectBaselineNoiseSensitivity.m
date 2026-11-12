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
learnrate = .05; % for adaptive mode 
if isEMU
    freqrng = [4, 9]; % Hz (Theta band)
else
    freqrng = [13, 30]; % Hz (beta band)
end

noiselevels = [0, .1*2.^(0:5)];

%% loop / parameter sweep 
errResults1 = cell(length(phTargets), width(dtaBL), length(noiselevels));
numResults1 = cell(size(errResults1));
trgResults1 = cell(size(errResults1));
errResults2 = errResults1; numResults2 = numResults1; trgResults2 = trgResults1;
% dim 1: target phase
% dim 2: channel 
% dim 3: noise level
% dim 4: const vs dynamic AR model 

for c = 1:1:width(dtaBL)

chtoplot = dtaBL.Properties.VariableNames{c}

for p = 1:length(phTargets)

phTarget = phTargets(p)

parfor l = 1:length(noiselevels)

noiselevel = noiselevels(l)

prog = ((c-1)*length(phTargets) + p-1 )*length(noiselevels) + l-1;
prog = prog/(length(phTargets)*width(dtaBL)*length(noiselevels));
disp([' ========== PROGRESS: ',num2str(round(100*prog)),'% ========== ']);

% run phase detection 
[phAll, phEst, frAll, frEst, trgTimeDyn, phStimDyn, ~, ~, durD, ...
    nCycleDyn, nMissDyn, nXtraDyn] = ...
    offline_PhaseDetect(dtaBL.(chtoplot)', Fs, [], dtaBL.Time', chtoplot, ...
    phTarget, freqrng, ARwin, ARord, predWin, -1, packetSize, ...
    0, false, false, false, false, false, true, noiselevel);
phErrDyn = radfix(phEst-phAll); frErrDyn = frEst - frAll;

pause(.001); 
drawnow;
pause(.001);

errResults1{p,c,l} = radfix(phStimDyn-phTarget);
trgResults1{p,c,l} = trgTimeDyn;
numResults1{p,c,l} = [nMissDyn, nXtraDyn, nCycleDyn];

[phAll, phEst, frAll, frEst, trgTimeDyn, phStimDyn, ~, ~, durD, ...
    nCycleDyn, nMissDyn, nXtraDyn] = ...
    offline_PhaseDetect(dtaBL.(chtoplot)', Fs, [], dtaBL.Time', chtoplot, ...
    phTarget, freqrng, ARwin, ARord, predWin, -1, packetSize, ...
    learnrate, false, false, false, false, false, true, noiselevel);
phErrDyn = radfix(phEst-phAll); frErrDyn = frEst - frAll;

pause(.001); 
drawnow;
pause(.001);

errResults2{p,c,l} = radfix(phStimDyn-phTarget);
trgResults2{p,c,l} = trgTimeDyn;
numResults2{p,c,l} = [nMissDyn, nXtraDyn, nCycleDyn];

end % noise levels
end % phase targets
end % channels

errResults = cat(4, errResults1, errResults2);
trgResults = cat(4, trgResults1, trgResults2);
numResults = cat(4, numResults1, numResults2);

%% aggregate/display all channel/target results 
mdlnames = {'Constant', 'Adaptive'};

errResultsAll = cell(size(errResults,4), size(errResults,3));
numResultsAll = cell(size(numResults,4), size(numResults,3));
for m = 1:size(errResults,4)
for r = 1:size(errResults,1)
    for h = 1:size(errResults,3)
        for c = 1:size(errResults,2)
            errResultsAll{m,h} = [errResultsAll{m,h}, errResults{r,c,h,m}];
            numResultsAll{m,h} = [numResultsAll{m,h}; numResults{r,c,h,m}];
        end
    end
end
for h = 1:size(numResults,3)
    numResultsAll{m,h} = sum(numResultsAll{m,h});
    numResultsAll{m,h} = [numResultsAll{m,h}(1:2), ...
        numResultsAll{m,h}(3)-sum(numResultsAll{m,h}(1:2))];
end
end

errResultsAll_avg = cellfun(@(ph) circ_mean(ph'), errResultsAll);
errResultsAll_std = cellfun(@(ph) circ_std(ph'),  errResultsAll);

fig1 = figure; 
for m = 1:size(errResultsAll,1)
subplot(1,size(errResultsAll,1),m)
for h = 1:size(errResultsAll,2)
    polarhistogram(errResultsAll{m,h}, 18); hold on; 
end
title(mdlnames{m})
subtitle('Phase Error (causal - offline)'); 
lgd = legend(string(1./noiselevels.^2)); lgd.Title.String = 'SNR';
end

fig2 = figure; 
for m = 1:size(errResultsAll,1)
ax(m) = subplot(size(errResultsAll,1),1,m);
bar(noiselevels, errResultsAll_avg(m,:)); 
hold on; grid on; 
errorbar(noiselevels,errResultsAll_avg(m,:), ...
    errResultsAll_std(m,:),errResultsAll_std(m,:), '.');
xlabel('Noise Levels'); ylabel('Phase Error (rad)'); 
legend('Circ Mean', '±1 Circ SD');
title(mdlnames{m})
end
linkaxes(ax, 'xy')

%% save 
thisfilever = getFileVersion(thisfilename);
[~,thisfilename] = fileparts(thisfilename);
svname = [fp,fn,'_',thisfilename,'_',thisfilever];

save(svname, 'errResults', 'numResults', 'trgResults', ...
    'phTargets', 'chselName', 'packetSize', 'ARwin', 'ARord', ...
    'learnrate', 'freqrng', 'noiselevels');
saveas(fig1, [svname,'_polhisto'], 'png');
saveas(fig2, [svname,'_rectbar'], 'png');