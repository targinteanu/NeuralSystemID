function [A, fig] = plotModelFit(EEGlist, EpocList, fitfun, plotchan)
% Fit a model to data and plot and return the results. 
% Inputs: 
%   EEGlist: array of EEGlab EEG structs corresponding to trials 
%   EpocList: Cell array of epocs to fit; same size as EEGlist, but cells
%             can have arrays of epocs of any length 
%   fitfun: function that takes in (train_data, test_data) as timetables
%           and outputs: 
%                   predicted train data [timetable] 
%                   predicted test data [timetable] 
%                   train error [struct with field pRMSE] 
%                   test error [ibid] 
%                   A matrix list [3D array] 
%   plotchan: name of channel to plot 

lnspc = {'r','b','m';
         '--',':','-';
         1.25,1.5,1}; 
lnidx = 1;
minV = inf; maxV = -inf;
A = [];
lgd = {'true'};
fig(1) = figure; hold on; 

for trl = 1:length(EEGlist)
    eeg = EEGlist(trl); 
    minV = min(minV, min(eeg.data(:)));
    maxV = max(maxV, max(eeg.data(:)));
    plot(eeg2timetable(eeg), plotchan, 'Color','k', 'LineWidth',2);
    curEpochs = EpocList{trl};
    curEpochs = [curEpochs, curEpochs(1)];
    for ep = 1:(length(curEpochs)-1)
        if lnidx > length(lnspc)
            lnidx = 1;
        end

        eeg = curEpochs(ep); eeg2 = curEpochs(ep+1);
        [trnPred, tstPred, trnE, tstE, Atrl] = ...
            fitfun(eeg2timetable(eeg), eeg2timetable(eeg2));
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

fig(2) = figure('Units','normalized', 'Position',[.05,.1,.9,.5]); 
subplot(121); heatmap(mean(A,3)); title('A mean'); 
subplot(122); heatmap(std(A,[],3)); title('A SD');

end