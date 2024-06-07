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
    epochT = 20; % s
    epoch_dt = 20; % s

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

plotchan = 'CZ'; 
lnspc = {'r','b','m';
         '--',':','-';
         1.25,1.5,1}; 
lnidx = 1;
minV = inf; maxV = -inf;
A = [];
lgd = {'true'};
figure; hold on; 

for trl = 1:length(curEEGlist)
    eeg = curEEGlist(trl); 
    minV = min(minV, min(eeg.data(:)));
    maxV = max(maxV, max(eeg.data(:)));
    plot(eeg2timetable(eeg), plotchan, 'Color','k', 'LineWidth',2);
    curEpochs = EpocList{trl};
    for ep = 1:(length(curEpochs)-1)
        if lnidx > length(lnspc)
            lnidx = 1;
        end

        eeg = curEpochs(ep); eeg2 = curEpochs(ep+1);
        [trnPred, tstPred, trnE, tstE, Atrl] = ...
            fitLTIauton(eeg2timetable(eeg), eeg2timetable(eeg2));
        A = cat(3,A,Atrl);

        plot(trnPred, plotchan, 'LineStyle',lnspc{2,lnidx}, ...
                                'Color',lnspc{1,lnidx}, ...
                                'LineWidth',lnspc{3,lnidx}); 
        plot(tstPred, plotchan, 'LineStyle',lnspc{2,lnidx}, ...
                                'Color',lnspc{1,lnidx}, ...
                                'LineWidth',lnspc{3,lnidx});

        lgd = [lgd, [num2str(100*(trnE.pRMSE),3),'% RMSE'], ...
                    [num2str(100*(tstE.pRMSE),3),'% RMSE']]; 

        lnidx = lnidx+1;
    end
end

grid on;
legend(lgd); 
Vrng = maxV - minV;
maxV = maxV + .5*Vrng; minV = minV - .5*Vrng; 
Vrng0 = diff(ylim);
if Vrng0 > Vrng
    ylim([minV,maxV]);
end

figure('Units','normalized', 'Position',[.05,.1,.9,.5]); 
subplot(121); heatmap(mean(A,3)); title('A mean'); 
subplot(122); heatmap(std(A,[],3)); title('A SD');