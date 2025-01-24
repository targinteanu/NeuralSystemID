%% Parkinson's Disease (PD) Project - train autonomous systems
% Train various state space systems to capture dynamics of cortical
% recordings in Parkinson's Disease subjects. 
% Autonomous only: does not include brain stimulation. 

%% load data file 
[fn,fp] = uigetfile('trainednet*.mat');
load(fullfile(fp,fn), 'bgNSS', 'dataTrain', 'dataTest');
disp([fp,' --- ',fn]);

%% validation params 
chdisp = [1; 10; 20]; chdisp = [chdisp; chdisp+width(dataTrain)/2];
kstep = .25; % s
kstep = ceil(kstep * dataTrain.Properties.SampleRate); % sample
Lval = 1000; % sample

dataTrainVal = dataTrain(1:Lval,:); dataTestVal = dataTest{1}(1:Lval,:);
H = height(chdisp);
figure; 
for p = 1:H
    subplot(H,2, 2*(p-1)+1);
    plottbl(dataTrainVal, chdisp(p), 'k',2);
    hold on; grid on;
    subplot(H,2, 2*(p-1)+2);
    plottbl(dataTestVal, chdisp(p), 'k', 2);
    hold on; grid on;
end

%% NSS (display only)
disp('Neural State Space - Training Validation')
bgNSStrain = myPredict(bgNSS, dataTrainVal, kstep, true);
disp('Neural State Space - Testing Validation')
bgNSStest = myPredict(bgNSS, dataTestVal, kstep, true);
for p = 1:H
    subplot(H,2, 2*(p-1)+1);
    plottbl(bgNSStrain, chdisp(p));
    hold on; grid on;
    subplot(H,2, 2*(p-1)+2);
    plottbl(bgNSStest, chdisp(p));
    hold on; grid on;
end

%% null system 
A = eye(width(dataTrain)); 
B = zeros(height(A),0); C = eye(size(A)); D = zeros(height(C),0);
sysNull = idss(ss(A,B,C,D, seconds(dataTrain.Properties.TimeStep)));
sysNull.StateName = dataTrain.Properties.VariableNames; 
sysNull.StateUnit = dataTrain.Properties.VariableUnits;
sysNull.OutputName = dataTrain.Properties.VariableNames; 
sysNull.OutputUnit = dataTrain.Properties.VariableUnits;

disp('Null State Space - Training Validation')
sysNulltrain = myPredict(sysNull, dataTrainVal, kstep, true);
disp('Null State Space - Testing Validation')
sysNulltest = myPredict(sysNull, dataTestVal, kstep, true);
for p = 1:H
    subplot(H,2, 2*(p-1)+1);
    plottbl(sysNulltrain, chdisp(p), ':k', 1.5);
    hold on; grid on;
    subplot(H,2, 2*(p-1)+2);
    plottbl(sysNulltest, chdisp(p), ':k', 1.5);
    hold on; grid on;
end

%% basic LTI system 
% full state = output
[~,~,~,~,A] = fitLTIauton(dataTrain);
B = zeros(height(A),0); C = eye(size(A)); D = zeros(height(C),0);
disp('Basic LTI - Training')
sysLTI = idss(ss(A,B,C,D, seconds(dataTrain.Properties.TimeStep)));
sysLTI.StateName = dataTrain.Properties.VariableNames; 
sysLTI.StateUnit = dataTrain.Properties.VariableUnits;
sysLTI.OutputName = dataTrain.Properties.VariableNames; 
sysLTI.OutputUnit = dataTrain.Properties.VariableUnits;

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

%% AR model
for p = 1:H
    disp(['AR - Training Channel ',num2str(chdisp(p))])
    arp = ar(dataTrain(:,chdisp(p)), 10, 'yw');
    disp('AR - Training Validation')
    xTrainPred = myPredict(arp, dataTrainVal(:,chdisp(p)), kstep, true);
    disp('AR - Testing Validation')
    xTestPred = myPredict(arp, dataTestVal(:,chdisp(p)), kstep, true);
    subplot(H,2, 2*(p-1)+1);
    plottbl(xTrainPred);
    hold on; grid on;
    subplot(H,2, 2*(p-1)+2);
    plottbl(xTestPred);
    hold on; grid on;
end

%% nontrivial LTI system 
N = 1;
LayerSize = [
    70;  % cortex 
    70;  % striatum 
    30;  % GP input 
    4;   % indirect input 
    1;   % STN input 
    6;   % STN output
    3;   % indirect output 
    10;  % GP output 
    250; % thalamus
    ];
%StateSize = N*sum(LayerSize);
StateSize = 200;
n4hzn = [1.5*StateSize, 7, 7];
disp('LTI - n4sid Training')
tic
bgLTI = n4sid(dataTrain, StateSize, ...
    n4sidOptions('Display','on', 'EstimateCovariance',false, ...
    'N4Weight','CVA', 'N4Horizon',n4hzn), ...
    'InputName',[],'OutputName',string(dataTrain.Properties.VariableNames));
toc
bgLTI.OutputName = dataTrain.Properties.VariableNames; 
bgLTI.OutputUnit = dataTrain.Properties.VariableUnits;

disp('LTI - Training Validation')
bgLTItrain = myPredict(bgLTI, dataTrainVal, kstep, false, true);
disp('LTI - Testing Validation')
bgLTItest = myPredict(bgLTI, dataTestVal, kstep, false, true);
for p = 1:H
    subplot(H,2, 2*(p-1)+1);
    plottbl(bgLTItrain, chdisp(p));
    hold on; grid on;
    subplot(H,2, 2*(p-1)+2);
    plottbl(bgLTItest, chdisp(p));
    hold on; grid on;
end

%% legend 
legend('true', 'bgNSS', 'null', 'sysLTI', 'AR', 'bgLTI')

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