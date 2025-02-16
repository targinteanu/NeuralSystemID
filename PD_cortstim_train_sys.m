%% Parkinson's Disease (PD) Project - train systems with input
% Train various state space systems to capture dynamics of cortical
% recordings in Parkinson's Disease subjects. 
% input = cortical brain stimulation. 

%% load data file 
[fn,fp] = uigetfile('sysLTI*.mat');
load(fullfile(fp,fn), 'dataStim');
disp([fp,' --- ',fn]);
[~,fn] = fileparts(fn);

%% divide stim into test-train 

fs = dataStim.Properties.SampleRate;
Nx = width(dataStim)-1;
OutputName = dataStim.Properties.VariableNames;
InputName = OutputName(end); OutputName = OutputName(1:(end-1));
OutputUnits = dataStim.Properties.VariableUnits(1:(end-1));

% reserve 4 min for training 
trainReserveDur = 4 * 60; % s
trainReserveN = ceil(trainReserveDur * fs);
dataTrain = dataStim(1:trainReserveN, :);
dataTest = dataStim((trainReserveN+1):end, :);

%% validation params 
chdisp = [1; 9; 18]; % chdisp = [chdisp; chdisp+width(dataTrain)/2];
%chdisp = [19; 38; 58];
kstep = .25; % s
kstep = ceil(kstep * dataTrain.Properties.SampleRate); % sample
Lval = 1000; % sample

dataTrainVal = dataTrain(1:Lval,:); dataTestVal = dataTest(1:Lval,:);
H = height(chdisp);
fig1 = figure('Units','normalized', 'Position',[.05,.1,.9,.8]); 
for p = 1:H
    ax(p,1) = subplot(H,2, 2*(p-1)+1);
    plottbl(dataTrainVal, chdisp(p), 'k',2);
    hold on; grid on;
    ax(p,2) = subplot(H,2, 2*(p-1)+2);
    plottbl(dataTestVal, chdisp(p), 'k', 2);
    hold on; grid on;
    linkaxes(ax(p,:), 'y');
end
linkaxes(ax(:,1), 'x'); linkaxes(ax(:,2), 'x');
subplot(H,2,1); title('Training'); subplot(H,2,2); title('Testing');

%% null system 
%{
A = eye(Nx); 
B = zeros(height(A),1); C = eye(size(A)); D = zeros(height(C),1);
sysNull = idss(ss(A,B,C,D, seconds(dataTrain.Properties.TimeStep)));
sysNull.StateName = OutputName; 
sysNull.StateUnit = OutputUnits;
sysNull.OutputName = OutputName; 
sysNull.OutputUnit = OutputUnits;

disp('Null State Space - Training Validation')
sysNulltrain = myPredict(sysNull, dataTrainVal(:,1:(end-1)), kstep, true);
disp('Null State Space - Testing Validation')
sysNulltest = myPredict(sysNull, dataTestVal(:,1:(end-1)), kstep, true);
for p = 1:H
    subplot(H,2, 2*(p-1)+1);
    plottbl(sysNulltrain, chdisp(p), ':k', 1.5);
    hold on; grid on;
    subplot(H,2, 2*(p-1)+2);
    plottbl(sysNulltest, chdisp(p), ':k', 1.5);
    hold on; grid on;
end
%}

%% basic LTI system 
%{
% full state = output
[~,~,~,~,A] = fitLTIauton(dataTrain);
B = zeros(height(A),0); C = eye(size(A)); D = zeros(height(C),0);
disp('Basic LTI - Training')
sysLTI = idss(ss(A,B,C,D, seconds(dataTrain.Properties.TimeStep)));
sysLTI.StateName = dataTrain.Properties.VariableNames; 
sysLTI.StateUnit = dataTrain.Properties.VariableUnits;
sysLTI.OutputName = dataTrain.Properties.VariableNames; 
sysLTI.OutputUnit = dataTrain.Properties.VariableUnits;

rat = sum([numel(sysLTI.A), numel(sysLTI.B), numel(sysLTI.C), numel(sysLTI.D), numel(sysLTI.K)]);
rat = numel(dataTrain)/rat; 
disp(['Training data is ',num2str(rat),' times parameter size'])

disp('Basic LTI - Training Validation')
sysLTItrain = myPredict(sysLTI, dataTrainVal, kstep, true);
disp('Basic LTI - Testing Validation')
sysLTItest = myPredict(sysLTI, dataTestVal, kstep, true);
for p = 1:H
    subplot(H,2, 2*(p-1)+1);
    plottbl(sysLTItrain, chdisp(p));
    hold on; grid on;
    subplot(H,2, 2*(p-1)+2);
    plottbl(sysLTItest, chdisp(p));
    hold on; grid on;
end
%}

%% AR model
%{
sysAR = [];
for p = 1:width(dataTrain)
    disp(['AR - Training Channel ',dataTrain.Properties.VariableNames{p}])
    ARp = ar(dataTrain(:,p), 10, 'yw');
    sysAR = [sysAR; ARp];
end

rat = height(dataTrain)/10; 
disp(['Training data is ',num2str(rat),' times parameter size'])

disp('AR - Training Validation')
xTrainPred = myPredict(sysAR, dataTrainVal, kstep, true);
disp('AR - Testing Validation')
xTestPred = myPredict(sysAR, dataTestVal, kstep, true);
for p = 1:H
    subplot(H,2, 2*(p-1)+1);
    plottbl(xTrainPred, chdisp(p));
    hold on; grid on;
    subplot(H,2, 2*(p-1)+2);
    plottbl(xTestPred, chdisp(p));
    hold on; grid on;
end
%}

%% nontrivial LTI system 
StateSize = 64;
n4hzn = [ceil(1.5*StateSize), 7, 7];
disp('LTI - n4sid Training')
tic
bgLTIstim = n4sid(dataTrain, StateSize, ...
    n4sidOptions('Display','on', 'EstimateCovariance',false, ...
    'N4Weight','CVA', 'N4Horizon',n4hzn), ...
    'InputName',InputName,'OutputName',OutputName);
toc
bgLTIstim.OutputName = OutputName; 
bgLTIstim.OutputUnit = OutputUnits;

rat = sum([numel(bgLTIstim.A), numel(bgLTIstim.B), numel(bgLTIstim.C), numel(bgLTIstim.D), numel(bgLTIstim.K)]);
rat = numel(dataTrain)/rat; 
disp(['Training data is ',num2str(rat),' times parameter size'])

plothelper(bgLTIstim, dataTrainVal, dataTestVal, kstep, chdisp);

%% hw - piecewise linear 
disp('HWpl - Training')
tic
bgHWpl = nlhw(dataTrain, bgLTIstim, 'idPiecewiseLinear', 'idPiecewiseLinear');
toc

plothelper(bgHWpl, dataTrainVal, dataTestVal, kstep, chdisp);

%% hw - sigmoid
disp('HWsg - Training')
tic
bgHWsg = nlhw(dataTrain, bgLTIstim, 'idSigmoidNetwork', 'idSigmoidNetwork');
toc

plothelper(bgHWsg, dataTrainVal, dataTestVal, kstep, chdisp);

%% hw - wavelet
disp('HWwl - Training')
tic
bgHWwl = nlhw(dataTrain, bgLTIstim, 'idWaveletNetwork', 'idWaveletNetwork');
toc

plothelper(bgHWwl, dataTrainVal, dataTestVal, kstep, chdisp);

%% hw - neural
disp('HWnn - Training')
tic
bgHWnn = nlhw(dataTrain, bgLTIstim, 'idNeuralNetwork', 'idNeuralNetwork');
toc

plothelper(bgHWnn, dataTrainVal, dataTestVal, kstep, chdisp);

%% legend 
legend('true', 'LTI', 'HWpl', 'HWsg', 'HWwl', 'HWnn')

%% saving 
svname = inputdlg('Save systems as:', 'File Save Name', 1, ...
    {[fn,'_andsyscortstimv']});
if ~isempty(svname)
    svname = svname{1};
    save(fullfile(fp,[svname,'.mat']), 'sysNull', 'sysLTI', 'sysAR', 'bgLTIstim', ...
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