%% Parkinson's Disease (PD) Project - train systems with input
% Train neural state space system to capture dynamics of cortical
% recordings in Parkinson's Disease subjects. 
% input = cortical brain stimulation. 

%% load data file 
[fn,fp] = uigetfile('sysLTI*cortstim*.mat'); % trained non-auton LTI, TF, and HW
load(fullfile(fp,fn), 'dataTest', 'dataTrain', 'bgLTI_NA'); 
disp([fp,' --- ',fn]);
[~,fn] = fileparts(fn);

%% validation params 
chdisp = [1; 9; 18]; % chdisp = [chdisp; chdisp+width(dataTrain)/2];
%chdisp = [19; 38; 58];
kstep = .25; % s
kstep = ceil(kstep * dataTrain.Properties.SampleRate); % sample
Lval = 1000; % sample

dataTrainVal = dataTrain(1:Lval,:); dataTestVal = dataTest(1:Lval,:);
C = height(chdisp);
fig1 = figure('Units','normalized', 'Position',[.05,.1,.9,.8]); 
for p = 1:C
    ax(p,1) = subplot(C,2, 2*(p-1)+1);
    plottbl(dataTrainVal, chdisp(p), 'k',2);
    hold on; grid on;
    ax(p,2) = subplot(C,2, 2*(p-1)+2);
    plottbl(dataTestVal, chdisp(p), 'k', 2);
    hold on; grid on;
    linkaxes(ax(p,:), 'y');
end
linkaxes(ax(:,1), 'x'); linkaxes(ax(:,2), 'x');
subplot(C,2,1); title('Training'); subplot(C,2,2); title('Testing');

%% build the model to be trained  

%ASize = size(bgLTI_NA.A); 
ASize = [64, 64];
OutputSize = width(dataTrain)-1;

% define a neural state space model 
bgNSS = idNeuralStateSpace(OutputSize, "NumInputs",1, "NumOutputs",OutputSize); 
bgNSS.StateNetwork = createMLPNetwork(bgNSS,"state", ...
    LayerSizes= ASize, ...
    WeightsInitializer="glorot", ...
    BiasInitializer="zeros", ...
    Activations='tanh');
bgNSS.StateName = dataTrain.Properties.VariableNames(1:(end-1)); 
bgNSS.StateUnit = dataTrain.Properties.VariableUnits(1:(end-1));
bgNSS.OutputName = dataTrain.Properties.VariableNames(1:(end-1)); 
bgNSS.OutputUnit = dataTrain.Properties.VariableUnits(1:(end-1));
bgNSS.InputName = dataTrain.Properties.VariableNames(end);

%% compare model size to data size 

numLearnables = [...
    bgNSS.StateNetwork.Learnables.Value; 
    bgNSS.OutputNetwork.Learnables.Value]; % cell array for each layer 
numLearnables = arrayfun(@(l) numel(l{:}), numLearnables); % total # at each layer
numLearnables = sum(numLearnables); % grand total # 

DataLearnableRatio = numel(dataTrain)/numLearnables;
disp(['Training data size is ',num2str(DataLearnableRatio),...
    ' times learnables size.']);

%% training 

%%{
% training options - ADAM
trnopts = nssTrainingOptions("adam");
trnopts.MaxEpochs = 1000; % default 100
trnopts.LearnRate = .001; % default .001
%trnopts.LearnRateSchedule = "piecewise"; % default "none"
trnopts.MiniBatchSize = 2048; % default 100
trnopts.LossFcn = "MeanSquaredError"; % default "MeanAbsoluteError"
%}

%{
% training options - SGDM
trnopts = nssTrainingOptions("sgdm");
trnopts.MaxEpochs = 1000; % default 100
trnopts.LearnRate = .0001; % default .01
%trnopts.LearnRateSchedule = "piecewise"; % default "none"
trnopts.MiniBatchSize = 4096; % default 1000
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
%{
bgNSS = nlssest(dataTrainEp, bgNSS, trnopts, ...
    UseLastExperimentForValidation = true);
%}
bgNSS = nlssest(dataTrain, bgNSS, trnopts);
toc

%% plot & legend
plothelper(bgNSS, dataTrainVal, dataTestVal, kstep, chdisp);
plothelper(bgLTI_NA, dataTrainVal, dataTestVal, kstep, chdisp);
legend('true', 'NSS', 'LTI');

%% saving 
svname = inputdlg('Save systems as:', 'File Save Name', 1, ...
    {[fn,'_andNNv']});
if ~isempty(svname)
    svname = svname{1};
    save(fullfile(fp,[svname,'.mat']), 'bgNSS', ...
        'dataTrain', 'dataTest', 'fn')
    saveas(fig1, fullfile(fp,svname),'fig'); 
    saveas(fig1, fullfile(fp,svname),'png'); 
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
    if ~isempty(TBL.Properties.VariableUnits)
        ylabel([TBL.Properties.VariableNames{v},' (',...
            TBL.Properties.VariableUnits{v},')']);
    else
        ylabel(TBL.Properties.VariableNames{v});
    end
    xlabel('time');
end

function [YPtrain, YPtest] = plothelper(sys, dataTrainVal, dataTestVal, kstep, chdisp)
disp(' - Training Validation')
YPtrain = predict(sys, dataTrainVal, kstep, predictOptions('InitialCondition','z'));
YPtrain.Time = YPtrain.Time + dataTrainVal.Time(1);
disp(' - Testing Validation')
YPtest = predict(sys, dataTestVal, kstep, predictOptions('InitialCondition','z'));
YPtest.Time = YPtest.Time + dataTestVal.Time(1);
H = height(chdisp);
for p = 1:H
    subplot(H,2, 2*(p-1)+1);
    plottbl(YPtrain, chdisp(p));
    hold on; grid on;
    subplot(H,2, 2*(p-1)+2);
    plottbl(YPtest, chdisp(p));
    hold on; grid on;
end
end