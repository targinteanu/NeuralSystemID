%% load EEG 
basedir = cd; 
cd('/Users/torenarginteanu/Desktop/Data_Chronic Pain'); 
[fn,fp] = uigetfile('.mat', 'Select Pre- or Post-processed EEG');
cd(basedir); 
load([fp,filesep,fn]);

%% start EEGlab 
    eeglabpath = '/Applications/MATLAB_R2021b.app/toolbox/eeglab2022.0';
    addpath(eeglabpath);
    eeglab

%% epoch 
[curEEGlist, EpocList] = epochBaseline(EEG_table,...
    'BaselineOpen','before experiment',1,1);

%% plot fit
plotchan = 'CZ'; 
A = plotModelFit(curEEGlist, EpocList, @(tsTbl, trTbl) fitLTIauton(tsTbl, trTbl), plotchan, [13 30]);

%% plot source-sink
chanlocs = curEEGlist(1).chanlocs;
[srcness,snkness] = SourceSink(A, chanlocs);