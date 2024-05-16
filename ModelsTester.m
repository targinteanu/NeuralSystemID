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

curEEGlist = EEG_table.BaselineOpen('before experiment'); 
curEEGlist = curEEGlist{:}; 

    % Determine the epoch duration and overlap: 
    epochT = .5; % s
    epoch_dt = .5; % s

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

%% plot
plotchan = 'Cz'; 
colrs = {'r','b','m'}; colridx = 1;
minV = inf; maxV = -inf;
figure; hold on; 
for trl = 1:length(curEEGlist)
    eeg = curEEGlist(trl);
    minV = min(minV, min(eeg.data(:)));
    maxV = max(maxV, max(eeg.data(:)));
    plot(eeg2timetable(eeg), plotchan, 'Color','k', 'LineWidth',1.5);
    curEpochs = EpocList{trl};
    for ep = 1:length(curEpochs)
        if colridx > length(colrs)
            colridx = 1;
        end
        eeg = curEpochs(ep); 
        [trnPred, tstPred] = fitLTIauton(eeg2timetable(eeg));
        plot(trnPred, plotchan, 'LineStyle','-', 'Color',colrs{colridx}); 
        plot(tstPred, plotchan, 'LineStyle',':', 'Color',colrs{colridx});
        colridx = colridx+1;
    end
end
grid on;
Vrng = maxV - minV;
maxV = maxV + .5*Vrng; minV = minV - .5*Vrng; 
ylim([minV,maxV]);