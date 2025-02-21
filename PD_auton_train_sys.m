%% Parkinson's Disease (PD) Project - train autonomous systems
% Train various state space systems to capture dynamics of cortical
% recordings in Parkinson's Disease subjects. 
% Autonomous only: does not include brain stimulation. 

%% load data file 
[fn,fp] = uigetfile('sysLTI*.mat');
load(fullfile(fp,fn), 'dataTrain', 'dataTest');
disp([fp,' --- ',fn]);
[~,fn] = fileparts(fn);

%% validation params 
chdisp = [1; 9; 18]; % chdisp = [chdisp; chdisp+width(dataTrain)/2];
%chdisp = [19; 38; 58];
kstep = .25; % s
kstep = ceil(kstep * dataTrain.Properties.SampleRate); % sample
Lval = 100; % sample

dataTrainVal = dataTrain(1:Lval,:); dataTestVal = dataTest{1}(1:Lval,:);
H = height(chdisp);
fig1 = figure('Units','normalized', 'Position',[.05,.05,.9,.9]); 
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

%% NSS (display only)
%{
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
%}

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

%% AR model
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

%% nontrivial LTI system 
StateSize = 128;
n4hzn = [ceil(1.5*StateSize), 7, 7];
disp('LTI - n4sid Training')
tic
bgLTI = n4sid(dataTrain, StateSize, ...
    n4sidOptions('Display','on', 'EstimateCovariance',false, ...
    'N4Weight','CVA', 'N4Horizon',n4hzn), ...
    'InputName',[],'OutputName',string(dataTrain.Properties.VariableNames));
toc
bgLTI.OutputName = dataTrain.Properties.VariableNames; 
bgLTI.OutputUnit = dataTrain.Properties.VariableUnits;

rat = sum([numel(bgLTI.A), numel(bgLTI.B), numel(bgLTI.C), numel(bgLTI.D), numel(bgLTI.K)]);
rat = numel(dataTrain)/rat; 
disp(['Training data is ',num2str(rat),' times parameter size'])

disp('LTI - Training Validation')
bgLTItrain = myPredict(bgLTI, dataTrainVal, kstep, true);
disp('LTI - Testing Validation')
bgLTItest = myPredict(bgLTI, dataTestVal, kstep, true);
for p = 1:H
    subplot(H,2, 2*(p-1)+1);
    plottbl(bgLTItrain, chdisp(p));
    hold on; grid on;
    subplot(H,2, 2*(p-1)+2);
    plottbl(bgLTItest, chdisp(p));
    hold on; grid on;
end

%% legend 
legend('true', 'null', 'sysLTI', 'AR', 'bgLTI')

%% saving 
svname = inputdlg('Save systems as:', 'File Save Name', 1, ...
    {[fn,'_andsysv']});
if ~isempty(svname)
    svname = svname{1};
    save(fullfile(fp,[svname,'.mat']), 'sysNull', 'sysLTI', 'sysAR', 'bgLTI', ...
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