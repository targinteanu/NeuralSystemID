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

curEEGlist = [EEG_table.PinPrick('before experiment'), EEG_table.PinPrick('after experiment')];
curEEGlist = [curEEGlist{:}]; 

    % Determine the epoch duration and overlap: 
    epochT = 1; % s
    epochT = .5*epochT;

EpocList = cell(size(curEEGlist));
for lstIdx = 1:length(curEEGlist)
    eeg = curEEGlist(lstIdx);
    if ~isempty(eeg)
        curEpochs3D = pop_epoch(eeg, {'11'}, [-epochT, epochT]);
        t = 1:curEpochs3D.trials; 
        curEpochs = repmat(curEpochs3D, size(t));
        for idx = 1:length(t)
            curEpoch = curEpochs(idx); 
            curEpoch.trials = 1;
            curEpoch.data = curEpoch.data(:,:,idx);
            curEpochs(idx) = curEpoch;
        end
        EpocList{lstIdx} = curEpochs;
    end
end

clear curEpoch curEpochs eeg idx lstIdx t 

%% plot fit
plotchan = 'CZ'; 
A = plotModelFit(curEEGlist, EpocList, @(tsTbl, trTbl) fitLTIauton(tsTbl, trTbl), plotchan);

%% plot source-sink
chanlocs = curEEGlist(1).chanlocs;
[srcness,snkness] = SourceSink(A, chanlocs);