%% Parkinson's Disease (PD) Project - evaluate systems with input
% not testing neural network
% input = cortical brain stimulation. 

%% load data file 
[fn,fp] = uigetfile('*andsyscortstim*.mat');
load(fullfile(fp,fn), 'bgLTI_A', 'bgLTI_NA2A', 'bgTF', 'bgLTI_NA', ...
        'dataTrain', 'dataTest');
disp([fp,' --- ',fn]);
[~,fn] = fileparts(fn);
sysName = {'bgLTI_A', 'bgLTI_NA2A', 'bgTF', 'bgLTI_NA'};
sysColr = {'b',       'c',          'm',    'r'};
sys = cellfun(@eval,sysName, 'UniformOutput',false);
fsNew = dataTrain.Properties.SampleRate;

%% select training/testing validation (subsets) 
%%{
Lval = 2000; % # samples 
disp([num2str(Lval/fsNew),' seconds selected for validation.'])
dataTest = dataTest(1:Lval,:); dataTrain = dataTrain(1:Lval,:);
%}

%% k-step ahead prediction 
% artifact duration is around 10ms, but we should verify different time
% scales. 
hznmkr = {'o', 'x', '+'}; % marker 
hznlwd = [2, 1.25, .75]; % line width
hzns = [.125, .5, 1.5]; % seconds  
hzns = ceil(hzns * fsNew); % # samples 

ypTrain = {}; ypTest = {};
for hzn = hzns
    disp(['Eval ',num2str(hzn),'-step']);
    %%{
    disp('Training Eval');
    ypTrain_ = cellfun(@(S) myPredict(S, dataTrain, hzn, true), sys, ...
        'UniformOutput',false);
    ypTrain = [ypTrain; ypTrain_];
    %}
    disp('Testing Eval');
    ypTest_ = cellfun(@(S) myPredict(S, dataTest, hzn, true), sys, ...
        'UniformOutput',false);
    ypTest = [ypTest; ypTest_];
end
clear ypTrain_ ypTest_ 

stimTest = dataTest(:,end); stimTrain = dataTrain(:,end); 
dataTest = dataTest(:,1:(end-1)); dataTrain = dataTrain(:,1:(end-1)); 

%% evaluation - assessment

alph = .05; % confidence level, uncorrected 
alph = alph/width(dataTrain); % Bonferroni correction 

pRMSE = @(yp, y) rmse(yp.Variables, y.Variables) ./ rms(y.Variables);

errsTest = nan(length(hzns), width(dataTrain), length(sys)); 
corsTest = errsTest; pcorsTest = corsTest;
errsTrain = errsTest; 
corsTrain = errsTrain; pcorsTrain = corsTrain;

for s = 1:length(sys)
    for k = 1:length(hzns)
        t1 = hzns(k) + 1;
        errsTrain(k,:,s) = pRMSE(ypTrain{k,s}(t1:end,:), dataTrain(t1:end,:));
        errsTest(k,:,s) = pRMSE(ypTest{k,s}(t1:end,:), dataTest(t1:end,:));
        for c = 1:width(dataTrain)
            %%{
            [corsTrain(k,c,s), pcorsTrain(k,c,s)] = ...
                corr(ypTrain{k,s}{t1:end,c}, dataTrain{t1:end,c}, 'Tail','right');
            %}
            [corsTest(k,c,s), pcorsTest(k,c,s)] = ...
                corr(ypTest{k,s}{t1:end,c}, dataTest{t1:end,c}, 'Tail','right');
        end
    end
end

figStem = figure('Units','normalized', 'Position',[.05,.05,.9,.9]);
K = length(hzns);
for k = 1:K
    kt = hzns(k)*1000/fsNew;
    kt = num2str(round(kt));

    subplot(K,3,3*(k-1)+1);
    for s = 1:length(sys)
        colr = sysColr{s};
        stem(100*errsTest(k,:,s), ['-',colr,'s'], 'LineWidth',1.5); hold on;
    end
    ylabel('% RMSE'); 
    title([kt,'ms prediction']);
    xticks(1:width(dataTrain)); xticklabels(dataTrain.Properties.VariableNames);
    xlabel('Channel Name')

    subplot(K,3,3*(k-1)+2);
    for s = 1:length(sys)
        colr = sysColr{s};
        stem(corsTest(k,:,s), ['-',colr,'s'], 'LineWidth',1.5); hold on;
    end
    ylabel('Pearsons \rho'); 
    title([kt,'ms prediction']);
    xticks(1:width(dataTrain)); xticklabels(dataTrain.Properties.VariableNames);
    xlabel('Channel Name')

    axp = subplot(K,3,3*(k-1)+3);
    for s = 1:length(sys)
        colr = sysColr{s};
        stem(pcorsTest(k,:,s), ['-',colr,'s'], 'LineWidth',1.5); hold on;
    end
    set(axp, 'YScale', 'log')
    xl = xlim(); plot(xl, alph*ones(size(xl)), 'k', 'LineWidth',1);
    ylabel('Correlation p-value'); 
    title([kt,'ms prediction']);
    xticks(1:width(dataTrain)); xticklabels(dataTrain.Properties.VariableNames);
    xlabel('Channel Name')

    %legend({'AR', 'LTI'}, "Location","best");

end

%% evaluation - visualization 

% best channel 
%corBest = min(corsTest, [], 3); % best model 
corBest = corsTest(:,:,end); % use bgLTI 
corRepr = abs(corBest - mean(corBest,2));
for k = 1:length(hzns)
    corBest(k,:) = assignrank(corBest(k,:)); % rank channel; higher = better 
    corRepr(k,:) = assignrank(corRepr(k,:)); % lower = more representative 
end
corBest = mean(corBest,1); % mean rank each horizon by channel 
corRepr = mean(corRepr,1);
[~,ch1] = max(corBest); [~,ch2] = min(corRepr);
% second-best 
selind = true(size(corBest)); selind(ch1) = false; selind(ch2) = false;
corBestSel = corBest(selind); corReprSel = corRepr(selind);
ch3 = max(corBestSel); ch4 = min(corReprSel);
ch3 = find(corBest == ch3); ch3 = ch3(1); 
ch4 = find(corRepr == ch4); ch4 = ch4(1); 
ch = [ch1; ch3; ch2; ch4];
%ch = unique(ch);
H = height(ch);

% plotting
figTime = figure('Units','normalized', 'Position',[.05,.05,.9,.9]);
for c = 1:H

    %%{
    % train - timeseries 
    ax(c,1) = subplot(H+1,4, 4*(c-1)+ 1);
    plottbl(dataTrain, ch(c), ':k', 4); grid on; hold on;
    for s = 1:width(ypTrain)
        for k = 1:length(hzns)
            plottbl(ypTrain{k,s}, ch(c), sysColr{s}, hznlwd(k)); 
        end
    end
    title('Train');

    % train - scatter 
    ax(c,2) = subplot(H+1,4, 4*(c-1)+ 2); hold on;
    for s = 1:width(ypTrain)
        colr = sysColr{s};
        for k = 1:length(hzns)
            mkr = hznmkr{k};
            plot(dataTrain{:,ch(c)}, ypTrain{k,s}{:,ch(c)}, [mkr,colr]);
        end
    end
    grid on; xlabel('true'); ylabel('pred');
    title('Train');
    %}

    % test - timeseries 
    ax(c,3) = subplot(H+1,4, 4*(c-1)+ 3);
    plottbl(dataTest, ch(c), ':k', 4); grid on; hold on;
    for s = 1:width(ypTest)
        for k = 1:length(hzns)
            plottbl(ypTest{k,s}, ch(c), sysColr{s}, hznlwd(k)); 
        end
    end
    title('Test');

    % test - scatter 
    ax(c,4) = subplot(H+1,4, 4*(c-1)+ 4); hold on;
    for s = 1:width(ypTest)
        colr = sysColr{s};
        for k = 1:length(hzns)
            mkr = hznmkr{k};
            plot(dataTest{:,ch(c)}, ypTest{k,s}{:,ch(c)}, [mkr,colr]);
        end
    end
    grid on; xlabel('true'); ylabel('pred');
    title('Test');

end

% stim (input) 
ax(H+1,1) = subplot(H+1,4, 4*H+ 1);
plottbl(stimTrain, 1, 'k', 1); grid on;
ax(H+1,3) = subplot(H+1,4, 4*H+ 3);
plottbl(stimTest, 1, 'k', 1); grid on;

linkaxes(ax(:,1), 'x'); linkaxes(ax(:,3), 'x');
linkaxes(ax(H+1,:), 'y');

%% sim visualization 
hzn = max(hzns);
ch = sort(ch);

nIR = 2; % # of impulse responses to show 

% training data; estimated initial condition 
dataTrainRT = retime([dataTrain, stimTrain], 'regular', 'nearest', 'TimeStep', seconds(bgLTIstim.Ts));
dataTrainRT.Time = dataTrainRT.Time - dataTrainRT.Time(1);
dataTestRT = retime([dataTest, stimTest], 'regular', 'nearest', 'TimeStep', seconds(bgLTIstim.Ts));
dataTestRT.Time = dataTestRT.Time - dataTestRT.Time(1);
ysLTI = sim(bgLTIstim, dataTrainRT(1:hzn, end), ...
    simOptions('InitialCondition', bgLTIstim.Report.Parameters.X0));

% training data; impulse response 
iStim = find(stimTrain.Variables > 0); %iStim = flipud(iStim);
mIR = floor(length(iStim)/nIR);
yIR_LTI_train = cell(2, nIR); nIdx = 1;
for idx = 1:mIR:length(iStim)
    idx_ = iStim(idx);
    it2 = ceil(.95*hzn) + idx_; it2 = min(it2, height(stimTrain));
    it1 = max(1, it2-hzn);
    x0 = pinv(bgLTIstim.C) * dataTrainRT{it1, 1:(end-1)}';
    yIR_LTI_train{1, nIdx} = sim(bgLTIstim, dataTrainRT(it1:it2, end), ...
        simOptions('InitialCondition', x0));
    yIR_LTI_train{2, nIdx} = dataTrainRT(it1:it2, :);
    yIR_LTI_train{1, nIdx}.Time = yIR_LTI_train{1, nIdx}.Time - yIR_LTI_train{1, nIdx}.Time(1);
    yIR_LTI_train{2, nIdx}.Time = yIR_LTI_train{2, nIdx}.Time - yIR_LTI_train{2, nIdx}.Time(1);
    nIdx = nIdx + 1;
end

% testing data; impulse response 
iStim = find(stimTest.Variables > 0); %iStim = flipud(iStim);
mIR = floor(length(iStim)/nIR);
yIR_LTI_test = cell(2, nIR); nIdx = 1;
for idx = 1:mIR:length(iStim)
    idx_ = iStim(idx);
    it2 = ceil(.95*hzn) + idx_; it2 = min(it2, height(stimTest));
    it1 = max(1, it2-hzn);
    x0 = pinv(bgLTIstim.C) * dataTestRT{it1, 1:(end-1)}';
    yIR_LTI_test{1, nIdx} = sim(bgLTIstim, dataTestRT(it1:it2, end), ...
        simOptions('InitialCondition', x0));
    yIR_LTI_test{2, nIdx} = dataTestRT(it1:it2, :);
    yIR_LTI_test{1, nIdx}.Time = yIR_LTI_test{1, nIdx}.Time - yIR_LTI_test{1, nIdx}.Time(1);
    yIR_LTI_test{2, nIdx}.Time = yIR_LTI_test{2, nIdx}.Time - yIR_LTI_test{2, nIdx}.Time(1);
    nIdx = nIdx + 1;
end

%{
% transform from log to uV 
dataTrainRT{:,1:(end-1)} = exp(dataTrainRT{:,1:(end-1)});
dataTestRT{:,1:(end-1)} = exp(dataTestRT{:,1:(end-1)});
ysLTI.Variables = exp(ysLTI.Variables);
for nIdx = 1:nIR
    yIR_LTI_train{1, nIdx}.Variables = exp(yIR_LTI_train{1, nIdx}.Variables);
    yIR_LTI_test{1, nIdx}.Variables = exp(yIR_LTI_test{1, nIdx}.Variables);
    yIR_LTI_train{2,nIdx}{:,1:(end-1)} = exp(yIR_LTI_train{2,nIdx}{:,1:(end-1)});
    yIR_LTI_test{2,nIdx}{:,1:(end-1)} = exp(yIR_LTI_test{2,nIdx}{:,1:(end-1)});
end
%}

% plot
Y = [{ysLTI; dataTrainRT(1:hzn,:)}, yIR_LTI_train, yIR_LTI_test];
W = width(Y);
figSim = figure('Units','normalized', 'Position',[.1,.1,.8,.8]);
for c = 1:H
    for n = 1:W
        ax2(c,n) = subplot(H+1,W, W*(c-1) +n);
        plottbl(Y{2,n}, ch(c), 'k', 3); hold on; grid on; 
        axis tight;
        plottbl(Y{1,n}, ch(c), sysColr{1}, 2);
        %legend('true', 'LTI');
        %set(ax2(c,n), 'YScale', 'log');
        set(ax2(c,n), "FontSize", 12);
    end
end
for n = 1:W
    ax2(H+1,n) = subplot(H+1,W, W*H +n); 
    stem(Y{2,n}.Time, Y{2,n}.stim, 'k');
    axis tight;
    set(ax2(H+1,n), "FontSize", 12);
    linkaxes(ax2(:,n), 'x');
end
%{
for c = 1:(H+1)
    linkaxes(ax2(c,:), 'y');
end
%}

%% saving 
svname = inputdlg('Save systems as:', 'File Save Name', 1, ...
    {[fn,'_eval']});
if ~isempty(svname)
    svname = svname{1};
    save(fullfile(fp,[svname,'.mat']), 'sys', 'hzns', 'sysName', ...
        'dataTrain', 'dataTest', 'ypTest', 'ypTrain', 'fn', ...
        'corsTest', 'pcorsTest', 'errsTest')
    saveas(figStem, fullfile(fp,[svname,'-stem']),'fig'); 
    saveas(figStem, fullfile(fp,[svname,'-stem']),'png'); 
    saveas(figTime, fullfile(fp,[svname,'-time']),'fig'); 
    saveas(figTime, fullfile(fp,[svname,'-time']),'png'); 
    saveas(figSim,  fullfile(fp,[svname,'-sim']),'fig'); 
    saveas(figSim,  fullfile(fp,[svname,'-sim']),'png'); 
end

%% helpers

function TBL = myPredictCell(sys, tblcell, k, b)
TBL = [];
for i = 1:length(tblcell)
    tbl = myPredict(sys, tblcell{i}, k, b);
    TBL = [TBL; tbl];
end
end

function r = assignrank(x)
[~,ord] = sort(x);
r = arrayfun(@(i) find(ord==i), 1:length(x));
end

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