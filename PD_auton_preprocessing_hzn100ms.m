%% Parkinson's Disease (PD) Project - preprocessing 
% Re-evaluate preprocessing using 100-ms horizon
% Autonomous only: does not include brain stimulation.

%% load the data 
[fn,fp] = uigetfile('sysLTIv*.mat');
load(fullfile(fp,fn), 'sysLTI', 'dataTrain', 'dataTest');
disp(fn);

[~,fn] = fileparts(fn);
svname = [fn,'_hzn100ms'];
fsNew = dataTrain.Properties.SampleRate;

%% visualize training results 
hzn = ceil(.1 * fsNew); % .25-second-ahead prediction horizon
Lval = 1000; % length of validation data 
ypTrain = myPredict(sysLTI, dataTrain(1:Lval,:), hzn, true);  
errsTrn = ...
    rmse(ypTrain.Variables, dataTrain(1:Lval,:).Variables) ./ ...
    rms(dataTrain(1:Lval,:).Variables);
errTrn = mean( errsTrn );
rhoTrn = arrayfun(@(c) corr(dataTrain{1:Lval,c}, ypTrain{:,c}), 1:width(dataTrain));
[~,ch1] = min(errsTrn); % best channel
[~,ch2] = min(abs(errTrn - errsTrn)); % most representative channel
fig1 = figure('Units','normalized', 'Position',[.05 .05 .9 .9]); 

subplot(3,2,1); 
plottbl(dataTrain(1:Lval,:), ch1); grid on; 
hold on; plottbl(ypTrain, ch1);
title('best channel');
subplot(3,2,3); 
plottbl(dataTrain(1:Lval,:), ch2); grid on; 
hold on; plottbl(ypTrain, ch2);
title('most representative channel');

subplot(3,2,2);
plot(dataTrain{1:Lval,ch1}, ypTrain{:,ch1}, '.');
xlabel('actual'); ylabel('predicted'); grid on;
title(dataTrain.Properties.VariableNames{ch1});
subplot(3,2,4);
plot(dataTrain{1:Lval,ch2}, ypTrain{:,ch2}, '.');
xlabel('actual'); ylabel('predicted'); grid on;
title(dataTrain.Properties.VariableNames{ch2});

subplot(3,1,3); stem(errsTrn); grid on; 
hold on; stem(rhoTrn);
xticks(1:width(ypTrain)); xticklabels(ypTrain.Properties.VariableNames);
xlabel('Channel Name'); ylabel('pRMSE / corr'); legend('RMS','corr');
title(['Train: Mean RMSE = ',num2str(100*errTrn),'%'])

if ~isempty(svname)
    sgtitle(svname);
    saveas(fig1, fullfile(fp,[svname,'_ResultsTrain']),'fig');
    saveas(fig1, fullfile(fp,[svname,'_ResultsTrain']),'png');
end

%% testing results 
dataTestCat = dataTest{1}(1:Lval,:);
ypTestCat = myPredict(sysLTI, dataTestCat, hzn);
fig2 = figure('Units','normalized', 'Position',[.05 .05 .9 .9]); 

subplot(3,2,1);
plottbl(dataTestCat, ch1); grid on;
hold on; plottbl(ypTestCat, ch1); 
title('best channel');
subplot(3,2,3); 
plottbl(dataTestCat, ch2); grid on;
hold on; plottbl(ypTestCat, ch2); 
title('most representative channel');

subplot(3,2,2);
plot(dataTestCat{:,ch1}, ypTestCat{:,ch1}, '.');
xlabel('actual'); ylabel('predicted'); grid on;
title(dataTestCat.Properties.VariableNames{ch1});
subplot(3,2,4);
plot(dataTestCat{:,ch2}, ypTestCat{:,ch2}, '.');
xlabel('actual'); ylabel('predicted'); grid on;
title(dataTestCat.Properties.VariableNames{ch2});

errsTst = ...
    rmse(ypTestCat.Variables, dataTestCat(1:Lval,:).Variables) ./ ...
    rms(dataTestCat(1:Lval,:).Variables);
errTst = mean( errsTst );
rhoTst = arrayfun(@(c) corr(dataTestCat{:,c}, ypTestCat{:,c}), 1:width(dataTestCat));
subplot(3,1,3); stem(errsTst); grid on; 
hold on; stem(rhoTst);
xticks(1:width(ypTestCat)); xticklabels(ypTestCat.Properties.VariableNames);
xlabel('Channel Name'); ylabel('pRMSE / corr'); legend('RMS','corr');
title(['Test: Mean RMSE = ',num2str(100*errTst),'%'])

if ~isempty(svname)
    sgtitle(svname);
    saveas(fig2, fullfile(fp,[svname,'_ResultsTest']),'fig');
    saveas(fig2, fullfile(fp,[svname,'_ResultsTest']),'png');
end

%% helpers 
function plottbl(TBL, v, lspc, lwid)
    if nargin < 4
        lwid = 1;
    end
    if nargin < 3
        lspc = '-';
    end
    if nargin < 2
        v = 1;
    end
    plot(TBL.Time, TBL{:,v}, lspc, 'LineWidth',lwid);
    ylabel([TBL.Properties.VariableNames{v},' (',...
        TBL.Properties.VariableUnits{v},')']);
    xlabel('time');
end