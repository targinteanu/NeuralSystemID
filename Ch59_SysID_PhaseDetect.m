%% load, define 

load("Ch59.mat", "NS2"); 
ind_bl_str = 4.5e5; ind_bl_end = 9.5e5; % no-stim baseline 
ind_rec_end = 2e6; % exclude dbs stim that is not annotated

Tr=NS2.Data(64,:);
Tr = Tr(1:ind_rec_end);
Tr_thr = Tr > 1e4;
nDelay = 3;
Tr_thr = [Tr_thr((nDelay+1):end), false(1,nDelay)]; % start nDelay samples earlier
Fs = NS2.MetaTags.SamplingFreq;

% noisy/ref channels are 3, 54, 57, 61
% include motor channels expected to have beta activity
chincl = {...
    'LS13','LS14','LS15','LS16','LS17',...
    'LS35','LS36','LS37','LS38','LS39','LS40','LS41',...
    'LS56','LS58','LS59','LS60','LS62'};
dta = ns2timetable(NS2); 
dta = dta(:,chincl); 
dta = dta(1:ind_rec_end,:); 

%% preprocess 

%Fs = dta.Properties.SampleRate;
hpf = designfilt('highpassiir', 'SampleRate',Fs, 'DesignMethod','butter', ...
    'StopbandFrequency',.5, 'PassbandFrequency',1.5, ...
    'StopbandAttenuation',60, 'PassbandRipple',1);
dta = FilterTimetable(@(d,x) filtfilt(d,x), hpf, dta);
dtaBL = dta(ind_bl_str:ind_bl_end,:);

packetSize = 10; 
phTarget = 0;
ARwin = 1000; ARlen = 50; predWin = 100;

errResults = cell(3, width(dtaBL), 2);
% dim 1: phase est err, freq est err, phase target err 
% dim 2: channel 
% dim 3: const vs dynamic AR model 

for c = 1:width(dtaBL)

chtoplot = dtaBL.Properties.VariableNames{c}

%% run phase detection 

% with constant AR model: 
[phAll, phEst, frAll, frEst, ~, phStimConst] = ...
    offline_PhaseDetect(dtaBL.(chtoplot)', Fs, [], dtaBL.Time', chtoplot, ...
    phTarget, [13,30], ARwin, ARlen, predWin, -1, packetSize, -1, [], true, true);
phErrConst = phEst-phAll; frErrConst = frEst - frAll;

pause(.001); 
drawnow;
pause(.001);

% with dynamic AR model: 
[phAll, phEst, frAll, frEst, ~, phStimDyn] = ...
    offline_PhaseDetect(dtaBL.(chtoplot)', Fs, [], dtaBL.Time', chtoplot, ...
    phTarget, [13,30], ARwin, ARlen, predWin, -1, packetSize, .0001, true, true, false);
phErrDyn = phEst-phAll; frErrDyn = frEst - frAll;

pause(.001); 
drawnow;
pause(.001);

%% compare constant vs dynamic AR results

figure; sgtitle(chtoplot);
subplot(1,2,1); 
polarhistogram(phErrConst); hold on; polarhistogram(phErrDyn);
title('Phase error (causal - offline)'); 
legend( ...
    ['Constant - RMSE ',num2str(rms(phErrConst))], ... 
    ['Dynamic - RMSE ',num2str(rms(phErrDyn))], ...
    'Location','northoutside');
subplot(1,2,2); 
histogram(frErrConst); hold on; grid on; histogram(frErrDyn);
xlabel('frequency (Hz)'); 
title('Freq. error (causal - offline)'); 
legend( ...
    ['Constant - RMSE ',num2str(rms(frErrConst))], ... 
    ['Dynamic - RMSE ',num2str(rms(frErrDyn))], ...
    'Location','northoutside');

pause(.001); 
drawnow;
pause(.001);

errResults{1,c,1} = phErrConst; errResults{1,c,2} = phErrDyn;
errResults{2,c,1} = frErrConst; errResults{2,c,2} = frErrDyn;
errResults{3,c,1} = phStimConst-phTarget; 
errResults{3,c,2} = phStimDyn-phTarget;

end

%% aggregate all channel results 

errResultsAll = cell(size(errResults,1), size(errResults,3));
for r = 1:size(errResults,1)
    for h = 1:size(errResults,3)
        for c = 1:size(errResults,2)
            errResultsAll{r,h} = [errResultsAll{r,h}, errResults{r,c,h}];
        end
    end
end

figure; 
subplot(1,3,1); 
polarhistogram(errResultsAll{1,1}); hold on; polarhistogram(errResultsAll{1,2});
title('Phase error (causal - offline)'); 
legend( ...
    ['Constant - RMSE ',num2str(rms(errResultsAll{1,1}))], ... 
    ['Dynamic - RMSE ',num2str(rms(errResultsAll{1,2}))], ...
    'Location','northoutside');
[~,p] = ttest2(errResultsAll{1,1}.^2, errResultsAll{1,2}.^2, 'tail', 'right');
subtitle(['p = ',num2str(p)])
subplot(1,3,2); 
histogram(errResultsAll{2,1}); hold on; grid on; histogram(errResultsAll{2,2});
xlabel('frequency (Hz)'); 
title('Freq. error (causal - offline)'); 
legend( ...
    ['Constant - RMSE ',num2str(rms(errResultsAll{2,1}))], ... 
    ['Dynamic - RMSE ',num2str(rms(errResultsAll{2,2}))], ...
    'Location','northoutside');
[~,p] = ttest2(errResultsAll{2,1}.^2, errResultsAll{2,2}.^2, 'tail', 'right');
subtitle(['p = ',num2str(p)])
subplot(1,3,3); 
polarhistogram(errResultsAll{3,1}); hold on; polarhistogram(errResultsAll{3,2});
title('Stim error (causal - offline)'); 
legend( ...
    ['Constant - RMSE ',num2str(rms(errResultsAll{3,1}))], ... 
    ['Dynamic - RMSE ',num2str(rms(errResultsAll{3,2}))], ...
    'Location','northoutside');
[~,p] = ttest2(errResultsAll{3,1}.^2, errResultsAll{3,2}.^2, 'tail', 'right');
subtitle(['p = ',num2str(p)])