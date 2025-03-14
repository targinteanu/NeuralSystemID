%% Parkinson's Disease (PD) Project - train systems with input
% Train various state space systems to capture dynamics of cortical
% recordings in Parkinson's Disease subjects. 
% input = cortical brain stimulation. 

%% load data file 
[fnA,fp] = uigetfile('sysLTI*.mat'); % trained auton sys
fn_ = find(fnA == '_'); fn_ = fn_(1); 
fn = [fnA(1:fn_-1),'.mat'];
load(fullfile(fp,fn), 'dataStim');
load(fullfile(fp,fnA), 'sysAR', 'sysLTI', 'bgLTI', 'sysNull');
disp([fp,' --- ',fn,' & ',fnA]);
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

%% nontrivial LTI system 
%StateSize = 128;
%StateSize = floor(StateSize/width(dataTrain)) * width(dataTrain); 
StateSize = height(bgLTI.A);
    % make state space multiple of output space
n4hzn = [ceil(1.5*StateSize), 7, 7];
disp('LTI - n4sid Training')
tic
bgLTI_NA = n4sid(dataTrain, StateSize, ...
    n4sidOptions('Display','on', 'EstimateCovariance',false, ...
    'N4Weight','CVA', 'N4Horizon',n4hzn), ...
    'InputName',InputName,'OutputName',OutputName);
toc
bgLTI_NA.OutputName = OutputName; 
bgLTI_NA.OutputUnit = OutputUnits;

rat = sum([numel(bgLTI_NA.A), numel(bgLTI_NA.B), numel(bgLTI_NA.C), numel(bgLTI_NA.D), numel(bgLTI_NA.K)]);
rat = numel(dataTrain)/rat; 
disp(['Training data is ',num2str(rat),' times parameter size'])

plothelper(bgLTI_NA, dataTrainVal, dataTestVal, kstep, chdisp);

%% refine 
%{
disp('LTI - refining with ssest')
tic
bgLTI_NA2NA = ssest(dataTrain, bgLTI_NA);
toc
plothelper(bgLTI_NA2NA, dataTrainVal, dataTestVal, kstep, chdisp);
%}

%% retrained 
bgLTI_A = idss( ss(bgLTI.A, zeros(size(bgLTI_NA.B)), bgLTI.C, zeros(size(bgLTI_NA.D)), ...
    (bgLTI.Ts) ) );
bgLTI_A.OutputName = OutputName; 
bgLTI_A.OutputUnit = OutputUnits;
bgLTI_A.InputName = bgLTI_NA.InputName;
%bgLTI_A.Report.Parameters.X0 = bgLTI.Report.Parameters.X0;

bgLTI_NA2A = bgLTI_NA; 
bgLTI_NA2A.B = 0.*bgLTI_NA2A.B; 
bgLTI_NA2A.D = 0.*bgLTI_NA2A.D; 

%{
disp('LTI - Retraining autonomous model with input-output data')
tic
bgLTI_A2NA = ssest(dataTrain, bgLTI_A, ...
    ssestOptions('Display','on', 'EstimateCovariance',false, ...
    'InitializeMethod','n4sid', 'InitialState','zero', ...
    'SearchMethod','gna', ...
    'N4Weight','CVA', 'N4Horizon',n4hzn));
toc
%}

plothelper(bgLTI_A, dataTrainVal, dataTestVal, kstep, chdisp);
plothelper(bgLTI_NA2A, dataTrainVal, dataTestVal, kstep, chdisp);
%plothelper(bgLTI_A2NA, dataTrainVal, dataTestVal, kstep, chdisp);

%% transfer function from auton model
disp('tf - estimating')
tic
bgTF = tfest(dataTrain, floor(StateSize/width(dataTrain)), ...
    tfestOptions('Display','on', 'EstimateCovariance',false, ...
    'InitialCondition','estimate', 'SearchMethod','gna'), ...
    'InputName',InputName,'OutputName',OutputName, ...
    'Ts',seconds(dataTrain.Properties.TimeStep) );
toc

plothelper(bgTF, dataTrainVal, dataTestVal, kstep, chdisp);

%%{
%% hw - piecewise linear 
disp('HWpl - Training')
tic
bgHWpl = nlhw(dataTrain, bgTF, 'idPiecewiseLinear', 'idPiecewiseLinear');
toc

plothelper(bgHWpl, dataTrainVal, dataTestVal, kstep, chdisp);

%% hw - sigmoid
disp('HWsg - Training')
tic
bgHWsg = nlhw(dataTrain, bgTF, 'idSigmoidNetwork', 'idSigmoidNetwork');
toc

plothelper(bgHWsg, dataTrainVal, dataTestVal, kstep, chdisp);

%% hw - wavelet
disp('HWwl - Training')
tic
bgHWwl = nlhw(dataTrain, bgTF, 'idWaveletNetwork', 'idWaveletNetwork');
toc

plothelper(bgHWwl, dataTrainVal, dataTestVal, kstep, chdisp);

%% hw - neural
disp('HWnn - Training')
tic
bgHWnn = nlhw(dataTrain, bgTF, 'idNeuralNetwork', 'idNeuralNetwork');
toc

plothelper(bgHWnn, dataTrainVal, dataTestVal, kstep, chdisp);

%% hw - deep neural
%{
disp('HWdn - Training')
tic
bgHWdn = nlhw(dataTrain, bgTF, ...
    idNeuralNetwork(NetworkType="dlnetwork"), idNeuralNetwork(NetworkType="dlnetwork"));
toc

plothelper(bgHWdn, dataTrainVal, dataTestVal, kstep, chdisp);
%}

%% legend 
legend('true', 'LTI-NA', 'LTI-A', 'LTI-NA2A', 'TF', ...
    'HWpl', 'HWsg', 'HWwl', 'HWnn')

%% saving 
svname = inputdlg('Save systems as:', 'File Save Name', 1, ...
    {[fn,'_andsyscortstimv']});
if ~isempty(svname)
    svname = svname{1};
    save(fullfile(fp,[svname,'.mat']), 'sysNull', 'sysLTI', 'sysAR', 'bgLTI_A', ...
        'bgLTI_NA', 'bgLTI_A', 'bgLTI_NA2A', 'bgTF', ...
        'bgHWpl', 'bgHWsg', 'bgHWwl', 'bgHWnn', ... 'bgHWdn', ...
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