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
ARwin = 5000; predWin = 100; % samples
ARord = 16; % # coefficients 
learnrate = .05; % for adaptive mode 
freqrng = [13, 30]; % Hz (beta band)

%% loop 
errResults = cell(length(phTargets), width(dtaBL), 2);
numResults = cell(size(errResults));
trgResults = cell(size(errResults));
% dim 1: target phase
% dim 2: channel 
% dim 3: const AR model vs const freq

durAR = []; durSin = [];

for c = 1:1:width(dtaBL)

chtoplot = dtaBL.Properties.VariableNames{c}

for p = 1:length(phTargets)

phTarget = phTargets(p)

prog = (c-1)*length(phTargets) + p;
prog = prog/(length(phTargets)*width(dtaBL));
disp([' ========== PROGRESS: ',num2str(round(100*prog)),'% ========== ']);

% run phase detection 

% with constant AR model: 
[phAll, phEst, frAll, frEst, trgTimeAR, phStimAR, ~, ~, durA, ...
    nCycleAR, nMissAR, nXtraAR] = ...
    offline_PhaseDetect(dtaBL.(chtoplot)', Fs, [], dtaBL.Time', chtoplot, ...
    phTarget, freqrng, ARwin, ARord, predWin, -1, packetSize, ...
    -1, [], false, false, false, false, true);
phErrAR = radfix(phEst-phAll); frErrAR = frEst - frAll;
durAR = [durAR, durA];

pause(.001); 
drawnow;
pause(.001);

% with const freq sinusoid: 
[phAll, phEst, frAll, frEst, trgTimeSin, phStimSin, ~, ~, durS, ...
    nCycleSin, nMissSin, nXtraSin] = ...
    offline_PhaseDetect(dtaBL.(chtoplot)', Fs, [], dtaBL.Time', chtoplot, ...
    phTarget, freqrng, ARwin, ARord, predWin, -1, packetSize, ...
    -1, [], false, false, false, false, false);
phErrSin = radfix(phEst-phAll); frErrSin = frEst - frAll;
durSin = [durSin, durS];

pause(.001); 
drawnow;
pause(.001);

errResults{p,c,1} = radfix(phStimAR-phTarget); 
errResults{p,c,2} = radfix(phStimSin-phTarget);
trgResults{p,c,1} = trgTimeAR; trgResults{p,c,2} = trgTimeSin;
numResults{p,c,1} = [nMissAR, nXtraAR, nCycleAR];
numResults{p,c,2} = [nMissSin, nXtraSin, nCycleSin];

end % phase targets
end % channels

%% aggregate/display all channel/target results 

durAR = durAR(~isnan(durAR)); durSin = durSin(~isnan(durSin));
fspc = 'took avg. %f, min %f, max %f seconds';
fspf = @(t) sprintf(fspc, mean(t), min(t), max(t));
disp(['Constant AR ',fspf(durAR)]);
disp(['Constant Freq Sin ',fspf(durSin)]);

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

fig1 = figure; 
polarhistogram(errResultsAll{1}, 18); hold on; polarhistogram(errResultsAll{2}, 18);
title('Stim error (causal - offline)'); 
legend( ...
    ['Constant AR - RMSE ',num2str(rms(errResultsAll{1}))], ... 
    ['Constant Sin Freq - RMSE ',num2str(rms(errResultsAll{2}))], ...
    'Location','northoutside');
[~,p] = ttest2(errResultsAll{1}.^2, errResultsAll{2}.^2, 'tail', 'right');
subtitle(['p = ',num2str(p)])

fig2 = figure; 
tiledlayout(1,2,'TileSpacing','compact');
ax1 = nexttile;
pie(ax1, numResultsAll{1});
title('Num. Stim. - AR');
ax2 = nexttile;
pie(ax2, numResultsAll{2});
title('Num. Stim. - Sin');
lgd = legend({'Missing', 'Extra', 'Correct'});
lgd.Layout.Tile = 'east';

%% save 
thisfilever = getFileVersion(thisfilename);
[~,thisfilename] = fileparts(thisfilename);
svname = [fp,fn,'_',thisfilename,'_',thisfilever];

save(svname, 'errResults', 'durAR', 'durSin', 'numResults', 'trgResults', ...
    'phTargets', 'chselName', 'packetSize', 'ARwin', 'ARord', 'learnrate', 'freqrng');
saveas(fig1, [svname,'_error'], 'png');
saveas(fig2, [svname,'_number'], 'png');