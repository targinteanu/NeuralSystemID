% comparing constant versus dynamic AR model 

% file selection
thisfilename = mfilename("fullpath");
[fn,fp] = uigetfile('*SegmentData*.mat', 'Choose Segmented Data File');
SegmentedDataFullfile = fullfile(fp,fn);
load(SegmentedDataFullfile);
[~,fn,fe] = fileparts(fn);

%% parameters 
packetSize = 20; % samples
phTargets = 2*pi*[0, .5, .3186, .6520, .7332, .1075]; % radians
ARwin = 1000; predWin = 500; % samples
ARord = 50; % # coefficients 
learnrate = .05; % for adaptive mode 
freqrng = [13, 30]; % Hz (beta band)

%% baseline data 

trngstr = @(Tbl) string(min(Tbl.Time))+" - "+string(max(Tbl.Time));

% data selection 
tblBaseline = tblsBaseline{1};
disp("Baseline: "+trngstr(tblBaseline))
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

%% fit baseline AR mdls 
BPF = buildFIRBPF(Fs,freqrng(1),freqrng(2));
ARmdls = cell(2,width(dtaBL));
for c = 1:width(dtaBL)
    disp(['Fitting baseline AR ',num2str(c),' of ',num2str(width(dtaBL))])
    xBLc = dtaBL{:,c};
    ARmdls{1,c} = ar(iddata(xBLc(1:ceil(3*ARwin)),[],1/Fs),ARord,'yw');
    xBLc = filtfilt(BPF,1,xBLc);
    ARmdls{2,c} = ar(iddata(xBLc(1:ceil(3*ARwin)),[],1/Fs),ARord,'yw');
end

%% non-baseline data 

tblsToTest = [tblsMisc(:,1); tblsSrl(:,1)];

tblsToTestDescs = cellfun(@(T) ...
    string(T.Properties.Description)+": "+trngstr(T), tblsToTest);
[selidx,selmade] = listdlg("ListString",tblsToTestDescs, ...
    "SelectionMode","multiple", "ListSize",[500,300], ...
    "PromptString","Select Test Conditions");
if selmade
    tblsToTest = tblsToTest(selidx);
end

tblsToTest2 = {};
for Ti = 1:length(tblsToTest)
    tblsToTest_Ti = tblsToTest{Ti};
    inan_Ti = isnan(tblsToTest_Ti.Variables);
    inan_Ti = sum(inan_Ti,2);
    if sum(inan_Ti)
        tblsToTest3 = {};
        for tj = find(inan_Ti)'
            tblsToTest_Tj = tblsToTest_Ti(1:(tj-1),:);
            tblsToTest_Ti = tblsToTest_Ti((ti+1):end,:);
            tblsToTest3 = [tblsToTest3; {tblsToTest_Tj}];
            inan_Ti = isnan(tblsToTest_Ti.Variables);
            inan_Ti = sum(inan_Ti,2);
        end
    else
        tblsToTest3 = {tblsToTest_Ti};
    end
    tblsToTest2 = [tblsToTest2; tblsToTest3];
end
tblsToTest2Descs = cellfun(@(T) ...
    string(T.Properties.Description)+": "+trngstr(T), tblsToTest2);

%% loop 
errResults = cell(length(phTargets), width(dtaBL), 2, length(tblsToTest2));
numResults = cell(size(errResults));
trgResults = cell(size(errResults));
% dim 1: target phase
% dim 2: channel 
% dim 3: const vs dynamic AR model 
% dim 4: table tested

durConst = []; durDyn = [];

for Ti = 1:length(tblsToTest2)

dta = tblsToTest2{Ti};
disp([' ***** ',num2str(Ti),' of ',num2str(length(tblsToTest2)),': '...
    dta.Properties.Description,' ***** '])

for c = 1:1:length(chselName)

chtoplot = chselName{c}

for p = 1:length(phTargets)

phTarget = phTargets(p)

prog = (c-1)*length(phTargets) + p;
prog = prog/(length(phTargets)*length(chselName));
disp([' ========== PROGRESS: ',num2str(round(100*prog)),'% ========== ']);

% run phase detection 

% with constant AR model: 
[phAll, phEst, frAll, frEst, trgTimeConst, phStimConst, ~, ~, durC, ...
    nCycleConst, nMissConst, nXtraConst] = ...
    offline_PhaseDetect(dta.(chtoplot)', Fs, [], dta.Time', chtoplot, ...
    phTarget, freqrng, ARwin, ARmdls(:,c), predWin, -1, packetSize, ...
    -1, [], false, false, false, false);
phErrConst = radfix(phEst-phAll); frErrConst = frEst - frAll;
durConst = [durConst, durC];

pause(.001); 
drawnow;
pause(.001);

% with dynamic AR model: 
[phAll, phEst, frAll, frEst, trgTimeDyn, phStimDyn, ~, ~, durD, ...
    nCycleDyn, nMissDyn, nXtraDyn] = ...
    offline_PhaseDetect(dta.(chtoplot)', Fs, [], dta.Time', chtoplot, ...
    phTarget, freqrng, ARwin, ARmdls(:,c), predWin, -1, packetSize, ...
    learnrate, false, false, false, false, false);
phErrDyn = radfix(phEst-phAll); frErrDyn = frEst - frAll;
durDyn = [durDyn, durD];

pause(.001); 
drawnow;
pause(.001);

errResults{p,c,1,Ti} = radfix(phStimConst-phTarget); 
errResults{p,c,2,Ti} = radfix(phStimDyn-phTarget);
trgResults{p,c,1,Ti} = trgTimeConst; trgResults{p,c,2} = trgTimeDyn;
numResults{p,c,1,Ti} = [nMissConst, nXtraConst, nCycleConst];
numResults{p,c,2,Ti} = [nMissDyn, nXtraDyn, nCycleDyn];

end % phase targets
end % channels
end % tables

%% aggregate/display all channel/target results 

durConst = durConst(~isnan(durConst)); durDyn = durDyn(~isnan(durDyn));
fspc = 'took avg. %f, min %f, max %f seconds';
fspf = @(t) sprintf(fspc, mean(t), min(t), max(t));
disp(['Constant ',fspf(durConst)]);
disp(['Dynamic ',fspf(durDyn)]);

errResultsAll = cell(1, size(errResults,3));
numResultsAll = cell(1, size(numResults,3));
for r = 1:size(errResults,1)
    for h = 1:size(errResults,3)
        for c = 1:size(errResults,2)
            for Ti = 1:size(errResults,4)
                errResultsAll{1,h} = [errResultsAll{1,h}, errResults{r,c,h,Ti}];
                numResultsAll{1,h} = [numResultsAll{1,h}; numResults{r,c,h,Ti}];
            end
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
    ['Constant - RMSE ',num2str(rms(errResultsAll{1}))], ... 
    ['Dynamic - RMSE ',num2str(rms(errResultsAll{2}))], ...
    'Location','northoutside');
[~,p] = ttest2(errResultsAll{1}.^2, errResultsAll{2}.^2, 'tail', 'right');
subtitle(['p = ',num2str(p)])

fig2 = figure; 
tiledlayout(1,2,'TileSpacing','compact');
ax1 = nexttile;
pie(ax1, numResultsAll{1});
title('Num. Stim. - Constant');
ax2 = nexttile;
pie(ax2, numResultsAll{2});
title('Num. Stim. - Dynamic');
lgd = legend({'Missing', 'Extra', 'Correct'});
lgd.Layout.Tile = 'east';

%% save 
thisfilever = getFileVersion(thisfilename);
[~,thisfilename] = fileparts(thisfilename);
svname = [fp,fn,'_',thisfilename,'_',thisfilever];

save(svname, 'errResults', 'durConst', 'durDyn', 'numResults', 'trgResults', ...
    'phTargets', 'chselName', 'packetSize', 'ARwin', 'ARord', 'learnrate', ...
    'freqrng', 'tblsToTest2Descs');
saveas(fig1, [svname,'_error'], 'png');
saveas(fig2, [svname,'_number'], 'png');