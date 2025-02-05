%% Parkinson's Disease (PD) Project - preprocessing 
% Test out preprocessing using a simple LTI system. 
% To be selected for further system training. 
% Fast evaluation using 250-ms horizon
% Autonomous only: does not include brain stimulation.

%% load the data 
[fn,fp] = uigetfile('*_DataTimeTables.mat');
load(fullfile(fp,fn));
disp(fn);
disp(DataTimeTables(:,[1,3,4])); 
dataBaseline = DataTimeTables{1,2};

%% check/visualize (some of) the data
SD = std(dataBaseline); isOut = isoutlier(dataBaseline, 'mean');
[~,ord] = sort(SD.Variables);
ch1 = [ord(1:5), ord((end-4):end)];

%% preprocess data 

% filter freq range 
loco = 13; hico = 30;
fsOrig = dataBaseline.Properties.SampleRate; 

% filtering bound rules 
minfac         = 3;    % this many (lo)cutoff-freq cycles in filter
min_filtorder  = 15;   % minimum filter length

% filter order 
if loco>0
    filtord = minfac*fix(fsOrig/loco);
elseif hico>0
    filtord = minfac*fix(fsOrig/hico);
end
if filtord < min_filtorder
    filtord = min_filtorder;
end

filtwts = fir1(filtord, [loco, hico]./(fsOrig/2));
filtfun = @(b,x) filtfilt(b,1,x); 
dataBaseline = FilterTimetable(filtfun,filtwts,dataBaseline);

%% inst freq 
%{
[~,dataFreq] = instPhaseFreqTblSmooth(dataBaseline, [loco hico]);
%dataFreq.Variables = dataFreq.Variables - 20;
%dataFreq.Variables = tanh((dataFreq.Variables-20)/10);
%dataFreq = instfreq(dataBaseline);
%dataFreq = retime(dataFreq, dataBaseline.Time, "spline");
dataBaseline = dataFreq;
%}

%% envelope/power
%%{
%dataBaseline.Variables = log(max(eps, envelope(dataBaseline.Variables)));
dataBaseline.Variables = envelope(dataBaseline.Variables);
for c = 1:width(dataBaseline)
    dataBaseline.Properties.VariableNames{c} = ...
        [dataBaseline.Properties.VariableNames{c},' envelope'];
    %{
    dataBaseline.Properties.VariableUnits{c} = ...
        ['log ',dataBaseline.Properties.VariableUnits{c}];
    %}
end
% power and freq
% dataBaseline = [dataBaseline, dataFreq];
%}

%% wavelet transform
%{
%N = 24; % vaid
N = 18; % beyl 
N = N/2;
dataWavelet1 = downsampleTimetable(dataBaseline, 2); dataWavelet2 = dataWavelet1;
for c = 1:width(dataBaseline)
    [y,x] = dwt(dataBaseline{:,c}, 'beyl'); % or try 'vaid'
    dataWavelet1{:,c} = y(N:end); dataWavelet2{:,c} = x(N:end);
    dataWavelet1.Properties.VariableNames{c} = ...
        [dataWavelet1.Properties.VariableNames{c},' Approx'];
    dataWavelet2.Properties.VariableNames{c} = ...
        [dataWavelet2.Properties.VariableNames{c},' Detail'];
end
%dataBaseline = [dataWavelet1, dataWavelet2]; 
dataBaseline = dataWavelet2;
%}

%% custom convolution 
wCustom = [1.5, -1, .6, -.33, .22, -.1, .05, 0, -.01, .05];
N = length(wCustom);
for c = 1:width(dataBaseline)
    x = conv(wCustom, dataBaseline{:,c});
    dataBaseline{:,c} = x(1:(end-N+1));
end

%% downsample, but ensure above nyquist rate 
fsNew = 2.1*hico;
fsRatio = floor(dataBaseline.Properties.SampleRate/fsNew); 
fsNew = dataBaseline.Properties.SampleRate / fsRatio; 
disp(['Resampling from ',num2str(fsOrig),' to ',num2str(fsNew)]);
dataBaseline = downsampleTimetable(dataBaseline, fsRatio);
fsNew = dataBaseline.Properties.SampleRate; 
disp(['Resampled to ',num2str(fsNew),' Hz']);

%% setup testing and training data 

% reserve 7 min for training 
trainReserveDur = 7 * 60; % s
trainReserveN = trainReserveDur * fsOrig; % samples at Fs of isOut 
% reserve the region with fewest total outliers in all channels
isOutSum = sum(isOut,2);
trainNoise = movmean(isOutSum, trainReserveN); 
trainNoise = trainNoise((ceil(trainReserveN/2)):(end-ceil(trainReserveN/2)+1));
[~,trainStartInd] = min(trainNoise); 
trainEndInd = min(height(isOutSum), trainStartInd + floor(trainReserveN));
trainNoise = isOutSum(trainStartInd:trainEndInd);
% convert to ind of new sampling rate 
trainStartInd = trainStartInd*fsNew/fsOrig;
trainStartInd = max(1, floor(trainStartInd));
trainReserveN = ceil(trainReserveN*fsNew/fsOrig);
trainEndInd = min(height(dataBaseline), trainStartInd + trainReserveN);
dataTrain = dataBaseline(trainStartInd:trainEndInd, :); 
% for some reason the TimeStep property sometimes gets messed up 
dataTrain = retime(dataTrain, 'regular', 'nearest', 'SampleRate', fsNew);
%{
dataTest = {...
    retime(dataBaseline(1:(trainStartInd-1), :), 'regular', 'nearest', 'SampleRate', fsNew); ...
    retime(dataBaseline((trainEndInd+1):end, :), 'regular', 'nearest', 'SampleRate', fsNew)};
%}
dataTest = {...
    dataBaseline(1:(trainStartInd-1), :); ...
    dataBaseline((trainEndInd+1):end, :)};
disp(string(dataTrain.Time(1))+" to "+string(dataTrain.Time(end))+" reserved for training.");

% try 20 x numLearnables = # training samples 
numLearnables = width(dataTrain)^2; 
DataLearnableRatio = numel(dataTrain)/numLearnables;
disp(['Training data size is ',num2str(DataLearnableRatio),...
    ' times learnables size.']);

%% training 

% run training 
tic
[~,~,~,~,A] = fitLTIauton(dataTrain);
B = zeros(height(A),0); C = eye(size(A)); D = zeros(height(C),0);
sysLTI = idss(ss(A,B,C,D, seconds(dataTrain.Properties.TimeStep)));
sysLTI.StateName = dataTrain.Properties.VariableNames; 
sysLTI.StateUnit = dataTrain.Properties.VariableUnits;
sysLTI.OutputName = dataTrain.Properties.VariableNames; 
sysLTI.OutputUnit = dataTrain.Properties.VariableUnits;
%sysLTI = d2c(sysLTI);
toc

% save training params 
svname = inputdlg('Save Training params as:', 'File Save Name', 1, {'sysLTIv'});
if ~isempty(svname)
    svname = svname{1};
    save(fullfile(fp,[svname,'.mat']), 'sysLTI', 'dataTrain', 'dataTest', 'fn')
end

%% visualize training results 
hzn = ceil(.25 * fsNew); % .25-second-ahead prediction horizon
Lval = 1000; % length of validation data 
ypTrain = myPredict(sysLTI, dataTrain(1:Lval,:), hzn, true);  
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

if ~isempty(svname)
    sgtitle(svname);
    saveas(fig1, fullfile(fp,[svname,'_ResultsTrain']),'fig');
    saveas(fig1, fullfile(fp,[svname,'_ResultsTrain']),'png');
end

%% testing results 
dataTestCat = dataTest{1}(1:Lval,:);
ypTestCat = myPredict(sysLTI, dataTestCat, hzn);
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