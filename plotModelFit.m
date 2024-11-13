function [A, trainEval, testEval, fig] = plotModelFit(...
    EEGlist, EpocList, fitfun, plotchan, fbnd, doEnvelope)
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
%   plotchan: name of channel to plot. If empty, no plot shown. 
%   fbnd: [ low cutoff , high cutoff ] frequency (Hz) to use for band-pass
%         filtering. If omitted or empty, the signal will not be filtered. 
%   doEnvelope: (true/false) take the envelope of the signal before fitting
%               the model, only if the signal is band-pass filtered. If
%               omitted, default is true if the signal is band-pass
%               filtered. 
% Outputs: 
%   A: state transition "A" matrix/matrices returned by fitfun. If there
%      are multiple epochs/trials with different A matrices, they will be
%      concatenated along the third dimension, i.e. A(:,:,trl) is the A
%      matrix from trial trl. 
%   trainEval: training evaluation (goodness-of-fit) statistics, including
%              fractional RMSE (pRMSE)
%   testEval: testing accuracy statistics in same format as above 
%   fig: figure handle(s) 

if nargin < 4
    plotchan = '';
end
if nargin < 5
    fbnd = [];
end
dofilt = ~isempty(fbnd);
if nargin < 6
    doEnvelope = dofilt;
end

if dofilt
    fs = EEGlist(1).srate;
    bpf = designfilt('bandpassiir', ...
        'SampleRate',fs, ...
        'HalfpowerFrequency1',fbnd(1), 'HalfpowerFrequency2',fbnd(2), ...
        'FilterOrder',20, ...
        'DesignMethod', 'butter');
end

lnspc = {'r','b','m';
         '--',':','-';
         1.25,1.5,1}; 
lnidx = 1;
minV = inf; maxV = -inf;
A = []; 
trainEval = []; testEval = [];
lgd = {'true'};
if isempty(plotchan)
    fig = [];
else
    fig(1) = figure; hold on; 
end

for trl = 1:length(EEGlist)
    eeg = EEGlist(trl); 
    minV = min(minV, min(eeg.data(:)));
    maxV = max(maxV, max(eeg.data(:)));
    eeg_ = eeg2timetable(eeg);
    if dofilt
        eeg_ = FilterTimetable(@(d,x)filtfilt(d,x), bpf, eeg_);
        if doEnvelope
            eeg_.Variables = envelope(eeg_.Variables);
        end
    end
    if ~isempty(plotchan)
        plot(eeg_, plotchan, 'Color','k', 'LineWidth',2);
    end
    curEpochs = EpocList{trl};
    curEpochs = [curEpochs, curEpochs(1)];
    for ep = 1:(length(curEpochs)-1)
        if lnidx > length(lnspc)
            lnidx = 1;
        end

        eeg = curEpochs(ep); eeg2 = curEpochs(ep+1);
        eeg_ = eeg2timetable(eeg);
        if dofilt
            eeg_ = FilterTimetable(@(d,x)filtfilt(d,x), bpf, eeg_);
            if doEnvelope
                eeg_.Variables = envelope(eeg_.Variables);
            end
        end
        [trnPred, tstPred, trnE, tstE, Atrl] = ...
            fitfun(eeg_, eeg2timetable(eeg2));
        A = cat(3,A,Atrl);
        trainEval = [trainEval, trnE]; testEval = [testEval, tstE];

        if ~isempty(plotchan)
            plot(trnPred, plotchan, 'LineStyle',lnspc{2,lnidx}, ...
                                    'Color',lnspc{1,lnidx}, ...
                                    'LineWidth',lnspc{3,lnidx}); 
            plot(tstPred, plotchan, 'LineStyle',lnspc{2,lnidx}, ...
                                    'Color',lnspc{1,lnidx}, ...
                                    'LineWidth',lnspc{3,lnidx});
        end

        lgd = [lgd, [num2str(100*(trnE.pRMSE),3),'% RMSE'], ...
                    [num2str(100*(tstE.pRMSE),3),'% RMSE']]; 

        lnidx = lnidx+1;
    end
end

if ~isempty(plotchan)
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

end