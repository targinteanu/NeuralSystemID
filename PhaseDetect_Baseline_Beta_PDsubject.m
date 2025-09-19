% file selection
thisfilename = mfilename("fullpath");
[fn,fp] = uigetfile('*SegmentData*.mat', 'Choose Segmented Data File');
SegmentedDataFullfile = fullfile(fp,fn);
load(SegmentedDataFullfile);
[~,fn,fe] = fileparts(fn);

% data selection 
tblBaseline = tblsBaseline{1};
chsel = contains(tblBaseline.Properties.VariableDescriptions, ...
    'PrG'); % Precentral Gyrus
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

%% parameters 
packetSize = 20; % samples
phTargets = 2*pi*[0, .5, .3186, .6520, .7332, .1075]; % radians
ARwin = 1000; predWin = 500; % samples
ARord = 50; % # coefficients 
learnrate = .05; % for adaptive mode 
freqrng = [13, 30]; % Hz (beta band)

%% loop 
errResults = cell(length(phTargets), width(dtaBL), 2);
% dim 1: target phase
% dim 2: channel 
% dim 3: const vs dynamic AR model 

durConst = []; durDyn = [];

for c = 1:1:width(dtaBL)

chtoplot = dtaBL.Properties.VariableNames{c}

for p = 1:length(phTargets)

phTarget = phTargets(p)

prog = (c-1)*length(phTargets) + p;
prog = prog/(length(phTargets)*width(dtaBL));
disp([' ========== PROGRESS: ',num2str(round(100*prog)),'% ========== ']);

% run phase detection 

% with constant AR model: 
[phAll, phEst, frAll, frEst, ~, phStimConst, ~, ~, durC] = ...
    offline_PhaseDetect(dtaBL.(chtoplot)', Fs, [], dtaBL.Time', chtoplot, ...
    phTarget, freqrng, ARwin, ARord, predWin, -1, packetSize, ...
    -1, [], false, false, false, false);
phErrConst = radfix(phEst-phAll); frErrConst = frEst - frAll;
durConst = [durConst, durC];

pause(.001); 
drawnow;
pause(.001);

% with dynamic AR model: 
[phAll, phEst, frAll, frEst, ~, phStimDyn, ~, ~, durD] = ...
    offline_PhaseDetect(dtaBL.(chtoplot)', Fs, [], dtaBL.Time', chtoplot, ...
    phTarget, freqrng, ARwin, ARord, predWin, -1, packetSize, ...
    learnrate, false, false, false, false, false);
phErrDyn = radfix(phEst-phAll); frErrDyn = frEst - frAll;
durDyn = [durDyn, durD];

pause(.001); 
drawnow;
pause(.001);

errResults{p,c,1} = phStimConst-phTarget; 
errResults{p,c,2} = phStimDyn-phTarget;

end % phase targets
end % channels

%% aggregate/display all channel/target results 

durConst = durConst(~isnan(durConst)); durDyn = durDyn(~isnan(durDyn));
fspc = 'took avg. %f, min %f, max %f seconds';
fspf = @(t) sprintf(fspc, mean(t), min(t), max(t));
disp(['Constant ',fspf(durConst)]);
disp(['Dynamic ',fspf(durDyn)]);

for r = 1:size(errResults,1)
    for h = 1:size(errResults,3)
        for c = 1:size(errResults,2)
            errResults{r,c,h} = radfix(errResults{r,c,h});
        end
    end
end

errResultsAll = cell(1, size(errResults,3));
for r = 1:size(errResults,1)
    for h = 1:size(errResults,3)
        for c = 1:size(errResults,2)
            errResultsAll{1,h} = [errResultsAll{1,h}, errResults{r,c,h}];
        end
    end
end

fig1 = figure; 
polarhistogram(errResultsAll{1}, 18); hold on; polarhistogram(errResultsAll{2}, 18);
title('Stim error (causal - offline)'); 
legend( ...
    ['Constant - RMSE ',num2str(rms(errResultsAll{1}))], ... 
    ['Dynamic - RMSE ',num2str(rms(errResultsAll{2}))], ...
    'Location','northoutside');
[~,p] = ttest2(errResultsAll{1}.^2, errResultsAll{2}.^2, 'tail', 'right');
subtitle(['p = ',num2str(p)])

%% save 
thisfilever = getFileVersion(thisfilename);
[~,thisfilename] = fileparts(thisfilename);
svname = [fp,fn,'_',thisfilename,'_',thisfilever];

save(svname, 'errResults', 'durConst', 'durDyn', ...
    'phTargets', 'chselName', 'packetSize', 'ARwin', 'ARord', 'learnrate', 'freqrng');
saveas(fig1, svname, 'png');