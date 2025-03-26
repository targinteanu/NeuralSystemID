%% Parkinson's Disease (PD) Project - train AR-LTI hybrid system
% Train a hybrid AR-LTI system to capture dynamics of cortical
% recordings in Parkinson's Disease subjects. 
% Autonomous only: does not include brain stimulation. 

%% load data file 
[fn,fp] = uigetfile('sysLTI*.mat');
load(fullfile(fp,fn), 'dataTrain', 'dataTest', 'sysLTI');
disp([fp,' --- ',fn]);
[~,fn] = fileparts(fn);
fsNew = dataTrain.Properties.SampleRate;

%% organize the data into taps 

numTaps = 10; % analog of AR model order 

% output, i.e. tap zero
dataTrainTapped = dataTrain(numTaps:end,:); 
dataTestTapped = cell(size(dataTest));
for trl = 1:length(dataTest)
    dataTestTapped{trl}  = dataTest{trl}(numTaps:end,:);
end

for n = 1:(numTaps-1)
    dataTrainTap = dataTrain( (numTaps-n):(end-n), :); 
    dataTrainTap.Time = dataTrainTapped.Time; 
    for c = 1:width(dataTrainTap)
        dataTrainTap.Properties.VariableNames{c} = ...
            [dataTrainTap.Properties.VariableNames{c},'_tap',num2str(n)];
    end
    dataTrainTapped = [dataTrainTapped, dataTrainTap]; 
    for trl = 1:length(dataTest)
        dataTestTap  = dataTest{trl}( (numTaps-n):(end-n), :);
        dataTestTap.Time  = dataTestTapped{trl}.Time;
        for c = 1:width(dataTestTap)
            dataTestTap.Properties.VariableNames{c} = ...
                [dataTestTap.Properties.VariableNames{c},'_tap',num2str(n)];
        end
        dataTestTapped{trl} = [dataTestTapped{trl}, dataTestTap];
    end
end

%% find the least-squares solution A matrix 

% check training size compared to learnables 
numLearnables = width(dataTrain)*width(dataTrainTapped);
DataLearnableRatio = numel(dataTrainTapped)/numLearnables;
disp(['Training data size is ',num2str(DataLearnableRatio),...
    ' times learnables size.']);

% solve the equation Y = CA*X, i.e. XT*CAT=YT
XT = dataTrainTapped.Variables; YT = dataTrain.Variables; YT = YT(numTaps:end,:);
XT = XT(1:(end-1),:); YT = YT(2:end,:);
CAT = XT\YT; CA = CAT';

% fill in unsolved rows for A 
A = zeros(width(CA)); 
C = zeros(width(dataTrain), height(A));
for m = 1:width(dataTrain) % each solved row 
    for n = 1:numTaps
        r = (m-1)*numTaps + n; % row 
        if n == 1
            % tap 0 = solved row 
            A(r,:) = CA(m,:); 
        else
            A(r,r-1) = 1;
        end
    end
    C(m, (m-1)*numTaps + 1) = 1;
end

%% define ARLTI system 
B = zeros(height(A),0); D = zeros(height(C),0);
sysARLTI = idss(ss(A,B,C,D, seconds(dataTrain.Properties.TimeStep)));
sysARLTI.StateName = dataTrainTapped.Properties.VariableNames; 
sysARLTI.StateUnit = dataTrainTapped.Properties.VariableUnits;
sysARLTI.OutputName = dataTrain.Properties.VariableNames; 
sysARLTI.OutputUnit = dataTrain.Properties.VariableUnits;

%% train AR model
sysAR = [];
for p = 1:width(dataTrain)
    disp(['AR - Training Channel ',dataTrain.Properties.VariableNames{p}])
    ARp = ar(dataTrain(:,p), numTaps, 'yw');
    sysAR = [sysAR; ARp];
end

rat = height(dataTrain)/numTaps; 
disp(['Training data is ',num2str(rat),' times parameter size'])

%% visualize training results 
hzn = ceil(.25 * fsNew); % .25-second-ahead prediction horizon
Lval = 10000; % length of validation data 
ypTrain = myPredict(sysARLTI, dataTrainTapped(1:Lval,:), hzn, true);  
ypTrain = predict(sysARLTI, dataTrainTapped(1:Lval,:), ...
    predictOptions('InitialCondition', dataTrainTapped{1,:}'));
ypTrain.Time = ypTrain.Time + dataTrainTapped.Time(1);
ypTrain = ypTrain((hzn+1):end,:); 
errsTrn = ...
    rmse(ypTrain.Variables, dataTrain((hzn+1):Lval,:).Variables) ./ ...
    rms(dataTrain((hzn+1):Lval,:).Variables);
errTrn = mean( errsTrn );
rhoTrn = arrayfun(@(c) corr(dataTrain{(hzn+1):Lval,c}, ypTrain{:,c}), 1:width(dataTrain));
[~,ch1] = min(errsTrn); % best channel
[~,ch2] = min(abs(errTrn - errsTrn)); % most representative channel
fig1 = figure('Units','normalized', 'Position',[.05 .05 .9 .9]); 

subplot(3,2,1); 
plottbl(dataTrain((hzn+1):Lval,:), ch1); grid on; 
hold on; plottbl(ypTrain, ch1);
title('best channel');
subplot(3,2,3); 
plottbl(dataTrain((hzn+1):Lval,:), ch2); grid on; 
hold on; plottbl(ypTrain, ch2);
title('most representative channel');

subplot(3,2,2);
plot(dataTrain{(hzn+1):Lval,ch1}, ypTrain{:,ch1}, '.');
xlabel('actual'); ylabel('predicted'); grid on;
title(dataTrain.Properties.VariableNames{ch1});
subplot(3,2,4);
plot(dataTrain{(hzn+1):Lval,ch2}, ypTrain{:,ch2}, '.');
xlabel('actual'); ylabel('predicted'); grid on;
title(dataTrain.Properties.VariableNames{ch2});

subplot(3,1,3); stem(errsTrn); grid on; 
hold on; stem(rhoTrn);
xticks(1:width(ypTrain)); xticklabels(ypTrain.Properties.VariableNames);
xlabel('Channel Name'); ylabel('pRMSE / corr'); legend('RMS','corr');
title(['Train: Mean RMSE = ',num2str(100*errTrn),'%'])

%% testing results 
dataTestCat = dataTest{1}(1:Lval,:);
ypTestCat = myPredict(sysARLTI, dataTestCat, hzn);
dataTestCat = dataTestCat((hzn+1):end,:); ypTestCat = ypTestCat((hzn+1):end,:);
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
    rmse(ypTestCat.Variables, dataTestCat.Variables) ./ ...
    rms(dataTestCat.Variables);
errTst = mean( errsTst );
rhoTst = arrayfun(@(c) corr(dataTestCat{:,c}, ypTestCat{:,c}), 1:width(dataTestCat));
subplot(3,1,3); stem(errsTst); grid on; 
hold on; stem(rhoTst);
xticks(1:width(ypTestCat)); xticklabels(ypTestCat.Properties.VariableNames);
xlabel('Channel Name'); ylabel('pRMSE / corr'); legend('RMS','corr');
title(['Test: Mean RMSE = ',num2str(100*errTst),'%'])

%% save 
svname = inputdlg('Save tapped system as:', 'File Save Name', 1, {'sysARLTIv'});
if ~isempty(svname)
    svname = svname{1};
    save(fullfile(fp,[svname,'.mat']), 'sysLTI', 'sysAR', 'sysARLTI', ...
        'dataTrain', 'dataTest', ...
        'dataTrainTapped', 'dataTestTapped', 'fn')
    figure(fig1); sgtitle(svname);
    saveas(fig1, fullfile(fp,[svname,'_ResultsTrain']),'fig');
    saveas(fig1, fullfile(fp,[svname,'_ResultsTrain']),'png');
    figure(fig2); sgtitle(svname); 
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