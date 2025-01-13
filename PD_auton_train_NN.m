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
minfac         = 1;    % this many (lo)cutoff-freq cycles in filter
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

% downsample, but ensure above nyquist rate 
fsNew = ceil(2.1*hico);
dataBaseline = resample(dataBaseline, fsNew, round(fsOrig));

%% build the model to be trained  

% define a deep neural network 
[~,bgDNN,numLearnables] = DLN_BasalGanglia_v2(width(dataBaseline), 1, true);
disp([num2str(numLearnables),' Learnables']);

% define a neural state space model 
bgNSS = idNeuralStateSpace(width(dataBaseline)); % autonomous; state = output
bgNSS.StateNetwork = bgDNN; 
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

% training options 
trnopts = nssTrainingOptions("adam");
trnopts.MaxEpochs = 1000;

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
hzn = ceil(.25 * fsNew); % .25-second-ahead prediction horizon
ypTrain = predict(bgNSS, dataTrain, hzn); 
ypTrain.Properties.VariableNames = dataTrain.Properties.VariableNames;
ypTrain.Time = ypTrain.Time + dataTrain.Time(1);
figure; myStackedPlot(dataTrain, chDisp(1)); 
hold on; myStackedPlot(ypTrain, chDisp(1)); 
errTrn = rmse(ypTrain.Variables, dataTrain.Variables, 'all');
title(['Total RMSE = ',num2str(errTrn)])

%% testing results 
ypTest = cellfun(@(T) predict(bgNSS, T, hzn), dataTest, 'UniformOutput', false);
for ep = 1:height(ypTest)
    ypTest{ep}.Time = ypTest{ep}.Time + dataTest{ep}.Time(1);
end
ypTestCat = [ypTest{1}; ypTest{2}]; 
dataTestCat = [dataTest{1}; dataTest{2}];
ypTestCat.Properties.VariableNames = dataTestCat.Properties.VariableNames;
figure; myStackedPlot(dataTestCat, chDisp(1));
hold on; myStackedPlot(ypTestCat, chDisp(1)); 
errTst = rmse(ypTestCat.Variables, dataTestCat.Variables, 'all');
title(['Total RMSE = ',num2str(errTst)])