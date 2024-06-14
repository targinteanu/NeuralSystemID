%% start EEGlab 
    eeglabpath = '/Applications/MATLAB_R2021b.app/toolbox/eeglab2022.0';
    addpath(eeglabpath);
    eeglab

%% load EEG 
basedir = cd; 
cd('/Users/torenarginteanu/Desktop/Data_Chronic Pain'); 
[fn,fp] = uigetfile('.mat', 'Select Pre- or Post-processed EEG');
cd(basedir); 
load([fp,filesep,fn]);

%% epoch 
[curEEGlist, EpocList] = epochStim(EEG_table,'PinPrick',1);

%% plot fit
plotchan = 'CZ'; 
A = plotModelFit(curEEGlist, EpocList, @(tsTbl, trTbl) fitLTIauton(tsTbl, trTbl), plotchan);

%% plot source-sink
chanlocs = curEEGlist(1).chanlocs;
[srcness,snkness,figs] = SourceSink(A, chanlocs);
figure(figs(2)); 
sgtitle(fn);