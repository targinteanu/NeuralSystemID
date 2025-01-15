%% load the data 
[fn,fp] = uigetfile('*_DataTimeTables.mat');
load(fullfile(fp,fn));
disp(fn);
disp(DataTimeTables(:,[1,3,4])); 
dataBaseline = DataTimeTables{1,2};

%% check/visualize (some of) the data
[BetaPower, SD, isOut, numOut] = tblChannelSummary(dataBaseline, [13, 30]);
[~,ord] = sort(SD.Variables);
chDisp = [ord(1:5), ord((end-4):end)];
figure; myStackedPlot(dataBaseline,chDisp,isOut);

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
[~,dataBaseline] = instPhaseFreqTbl(dataBaseline);
%dataBaseline.Variables = envelope(dataBaseline.Variables);

% downsample, but ensure above nyquist rate 
fsNew = ceil(2.1*hico);
dataBaseline = resample(dataBaseline, fsNew, round(fsOrig));

%% build the model to be trained  

% define a deep neural network 
[~,bgDNN,numLearnables] = DLN_BasalGanglia_v2(width(dataBaseline), 1, true);
disp([num2str(numLearnables),' Learnables']);

%% define a neural state space model 
bgNSS = idNeuralStateSpace(width(dataBaseline)); % autonomous; state = output
%bgNSS.StateNetwork = bgDNN; 
bgNSS.StateNetwork = createMLPNetwork(bgNSS,"state", ...
    LayerSizes=[256 256], ...
    WeightsInitializer="glorot", ...
    BiasInitializer="zeros", ...
    Activations='tanh');
bgNSS.StateName = dataBaseline.Properties.VariableNames; 
bgNSS.StateUnit = dataBaseline.Properties.VariableUnits;
bgNSS.OutputName = dataBaseline.Properties.VariableNames; 
bgNSS.OutputUnit = dataBaseline.Properties.VariableUnits;

%% setup testing and training data 

% reserve 7 min for training 
trainReserveDur = 7 * 60; % s
trainReserveN = trainReserveDur * fsOrig; % samples at Fs of isOut 
% reserve the region with fewest total outliers in all channels
isOutSum = sum(isOut,2);
trainNoise = movmean(isOutSum, trainReserveN); 
trainNoise = trainNoise((ceil(trainReserveN/2)):(end-ceil(trainReserveN/2)+1));
[~,trainStartInd] = min(trainNoise); 
trainEndInd = min(height(isOutSum), trainStartInd + trainReserveN);
trainNoise = isOutSum(trainStartInd:trainEndInd);
% convert to ind of new sampling rate 
trainStartInd = trainStartInd*fsNew/fsOrig;
trainStartInd = max(1, floor(trainStartInd));
trainReserveN = ceil(trainReserveN*fsNew/fsOrig);
trainEndInd = min(height(dataBaseline), trainStartInd + trainReserveN);
dataTrain = dataBaseline(trainStartInd:trainEndInd, :); 
% for some reason the TimeStep property sometimes gets messed up 
dataTrain = retime(dataTrain, 'regular', 'nearest', 'SampleRate', fsNew);
dataTest = {...
    retime(dataBaseline(1:(trainStartInd-1), :), 'regular', 'nearest', 'SampleRate', fsNew); ...
    retime(dataBaseline((trainEndInd+1):end, :), 'regular', 'nearest', 'SampleRate', fsNew)};
disp(string(dataTrain.Time(1))+" to "+string(dataTrain.Time(end))+" reserved for training.");

% put training data into 1s non-overlapping epochs 
epochdur = 1; % s 
epochNorig = round(epochdur/fsNew); % length at orig sample rate 
epochN = round(epochdur*fsNew); % length at new sample rate
numEpochTrain = floor(height(dataTrain)/epochN);
epochNoise = zeros(numEpochTrain, 1);
dataTrainEp = cell(numEpochTrain, 1);
for ep = 1:numEpochTrain
    epInd1 = (ep-1)*epochNorig + 1; epInd2 = epInd1 + epochNorig - 1;
    epochNoise(ep) = sum(trainNoise(epInd1:epInd2));
    epInd1 = (ep-1)*epochN + 1;     epInd2 = epInd1 + epochN - 1;
    epTT = dataTrain(epInd1:epInd2, :);
    epTT.Time = epTT.Time - epTT.Time(1); % all start at "0" 
    epTT = retime(epTT, 'regular', 'nearest', 'SampleRate', fsNew);
    dataTrainEp{ep} = epTT;
end
[~,epInd] = sort(epochNoise);

% try 20 x numLearnables = # training samples 
numelTrain = 20 * numLearnables; 
nTrain = ceil(numelTrain / width(dataTrain)); % # of time samples 
NTrain = ceil(nTrain / epochN); % # of epochs 
NTrain = NTrain+1; % one to be used for validation 
epInd = epInd(1:NTrain); dataTrainEp = dataTrainEp(epInd);

%% training 

%%{
% training options - ADAM
trnopts = nssTrainingOptions("adam");
trnopts.MaxEpochs = 1000; % default 100
%trnopts.LearnRate = .01; % default .001
%trnopts.LearnRateSchedule = "piecewise"; % default "none"
trnopts.MiniBatchSize = 4096; % default 100
trnopts.LossFcn = "MeanSquaredError"; % default "MeanAbsoluteError"
%}

%{
% training options - SGDM
trnopts = nssTrainingOptions("sgdm");
trnopts.MaxEpochs = 1000; % default 100
trnopts.LearnRate = .001; % default .01
%trnopts.LearnRateSchedule = "piecewise"; % default "none"
trnopts.MiniBatchSize = 2048; % default 1000
trnopts.LossFcn = "MeanSquaredError"; % default "MeanAbsoluteError"
%}

%{
% training options - LBFGS; only since MATLAB 2024b
trnopts = nssTrainingOptions("lbfgs");
trnopts.MaxIterations = 1000; % default 100
%trnopts.LossFcn = "MeanSquaredError"; % default "MeanAbsoluteError"
%}

% run training 
tic
bgNSS = nlssest(dataTrainEp, bgNSS, trnopts, ...
    UseLastExperimentForValidation = true);
toc

% save trained net 
svname = inputdlg('Save Trained NSS as:', 'File Save Name', 1, {'trainednetv'});
svname = [svname{1},'.mat'];
save(fullfile(fp,svname), 'bgNSS', 'dataTrainEp', 'dataTrain', 'dataTest', 'fn')

% visualize training results 
%hzn = ceil(.25 * fsNew); % .25-second-ahead prediction horizon
Lval = 1000; % length of validation data 
hzn = .25; % .25-second-ahead prediction horizon
ypTrain = predict(bgNSS, dataTrain(1:Lval,:), hzn); 
ypTrain.Properties.VariableNames = dataTrain.Properties.VariableNames;
ypTrain.Properties.VariableUnits = dataTrain.Properties.VariableUnits;
ypTrain.Time = ypTrain.Time + dataTrain.Time(1);
figure; myStackedPlot(dataTrain(1:Lval,:), chDisp(1)); 
hold on; myStackedPlot(ypTrain, chDisp(1)); 
errTrn = mean(...
    rmse(ypTrain.Variables, dataTrain(1:Lval,:).Variables) ./ ...
    rms(dataTrain(1:Lval,:).Variables) );
title(['Mean RMSE = ',num2str(100*errTrn),'%'])

%% testing results 
%{
ypTest = cellfun(@(T) predict(bgNSS, T, hzn), dataTest, 'UniformOutput', false);
for ep = 1:height(ypTest)
    ypTest{ep}.Time = ypTest{ep}.Time + dataTest{ep}.Time(1);
end
ypTestCat = [ypTest{1}; ypTest{2}]; 
dataTestCat = [dataTest{1}; dataTest{2}];
%}
dataTestCat = dataTest{1}(1:Lval,:);
ypTestCat = predict(bgNSS, dataTestCat, hzn);
ypTestCat.Time = ypTestCat.Time + dataTestCat.Time(1);
ypTestCat.Properties.VariableNames = dataTestCat.Properties.VariableNames;
ypTestCat.Properties.VariableUnits = dataTestCat.Properties.VariableUnits;
figure; myStackedPlot(dataTestCat, chDisp(1));
hold on; myStackedPlot(ypTestCat, chDisp(1)); 
errTst = mean(...
    rmse(ypTestCat.Variables, dataTestCat(1:Lval,:).Variables) ./ ...
    rms(dataTestCat(1:Lval,:).Variables) );
title(['Mean RMSE = ',num2str(100*errTst),'%'])

%% confirm train/test with simulation 

simStart = 1; % start index of sim 
simOptsTrn = simOptions(...
    'InitialCondition', dataTrain{simStart,:}', ...
    'OutputTimes', seconds(dataTrain.Time(simStart:end)-dataTrain.Time(simStart)));
ysTrain = sim(bgNSS, [], simOptsTrn);
ysTrain = array2timetable(ysTrain, "RowTimes",dataTrain.Time(simStart:end), ...
    "VariableNames",dataTrain.Properties.VariableNames);
ysTrain.Properties.VariableUnits = dataTrain.Properties.VariableUnits;
figure; myStackedPlot(dataTrain, chDisp(1)); 
hold on; myStackedPlot(ysTrain, chDisp(1)); 
hold on; myStackedPlot(ypTrain, chDisp(1));