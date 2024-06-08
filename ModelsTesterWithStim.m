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

curEEGlist = [EEG_table.PinPrick('before experiment'), EEG_table.PinPrick('after experiment')];
curEEGlist = [curEEGlist{:}]; 

    % Determine the epoch duration and overlap: 
    epochT = 1; % s
    epochT = .5*epochT;

EpocList = cell(size(curEEGlist));
for lstIdx = 1:length(curEEGlist)
    eeg = curEEGlist(lstIdx);
    if ~isempty(eeg)
        t = (eeg.xmin):epoch_dt:((eeg.xmax)-epochT);
        curEpochs = repmat(eeg, size(t));
        for idx = 1:length(t)
            if ~mod(idx/length(t), .05)
                disp(['Epoch ',num2str(idx),' of ',num2str(length(t)),...
                    ' (',num2str(100*idx/length(t),3),'%)'])
            end
            curEpoch = pop_select(eeg, 'time', t(idx)+[0,epochT]);
            curEpoch.xmin = curEpoch.xmin + t(idx);
            curEpoch.xmax = curEpoch.xmax + t(idx);
            curEpoch.times = curEpoch.times + t(idx)*1000;
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