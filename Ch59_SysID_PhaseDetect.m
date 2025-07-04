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

chtoplot = 'LS59';

%% preprocess 

%Fs = dta.Properties.SampleRate;
hpf = designfilt('highpassiir', 'SampleRate',Fs, 'DesignMethod','butter', ...
    'StopbandFrequency',.5, 'PassbandFrequency',1.5, ...
    'StopbandAttenuation',60, 'PassbandRipple',1);
dta = FilterTimetable(@(d,x) filtfilt(d,x), hpf, dta);
dtaBL = dta(ind_bl_str:ind_bl_end,:);

%% run phase detection 

packetSize = 10; 
phTarget = 0;
ARwin = 1000; ARlen = 10; predWin = 100;

% with constant AR model: 
[phAll, phEst, frAll, frEst] = ...
    offline_PhaseDetect(dtaBL.(chtoplot)', Fs, [], dtaBL.Time', chtoplot, ...
    phTarget, [13,30], ARwin, ARlen, predWin, -1, packetSize, -1, [], true);
phErrConst = phEst-phAll; frErrConst = frEst - frAll;

pause(.001); 
drawnow;
pause(.001);

% with dynamic AR model: 
[phAll, phEst, frAll, frEst] = ...
    offline_PhaseDetect(dtaBL.(chtoplot)', Fs, [], dtaBL.Time', chtoplot, ...
    phTarget, [13,30], ARwin, ARlen, predWin, -1, packetSize, .001, true, true);
phErrDyn = phEst-phAll; frErrDyn = frEst - frAll;

%% compare constant vs dynamic AR results
figure; 
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