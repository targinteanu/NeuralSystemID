%% Parkinson's Disease (PD) Project - evaluate autonomous systems
% Autonomous only: does not include brain stimulation. 
% not testing neural network
% simplified for speed and plot clarity. Removed some systems and no
% training validation. 

%% load data file 
[fn,fp] = uigetfile('*andsys*.mat');
load(fullfile(fp,fn), 'sysNull', 'sysLTI', 'sysAR', 'bgLTI', ...
    'dataTrain', 'dataTest');
disp([fp,' --- ',fn]);
[~,fn] = fileparts(fn);
sysName = {'sysAR', 'bgLTI'};
sysColr = {'b',     'r'};
sys = cellfun(@eval,sysName, 'UniformOutput',false);
fsNew = dataTrain.Properties.SampleRate;

[~,dataTestSel] = max(cellfun(@height, dataTest));
dataTest = dataTest{dataTestSel};

%% select training/testing validation (subsets) 
%{
Lval = 1000; % # samples 
disp([num2str(Lval/fsNew),' seconds selected for validation.'])
dataTest = dataTest(1:Lval,:); dataTrain = dataTrain(1:Lval,:);
%}

%% k-step ahead prediction 
% artifact duration is around 10ms, but we should verify different time
% scales. 
hznmkr = {'o', 'x'}; % marker 
hznlwd = [2, 1.25]; % line width
hzns = [.5, 1.5]; % seconds 
hzns = ceil(hzns * fsNew); % # samples 

ypTrain = {}; ypTest = {};
for hzn = hzns
    disp(['Eval ',num2str(hzn),'-step']);
    %{
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

%% evaluation - assessment

pRMSE = @(yp, y) rmse(yp.Variables, y.Variables) ./ rms(y.Variables);

errsTest = nan(length(hzns), width(dataTrain), length(sys)); 
corsTest = errsTest; pcorsTest = corsTest;
%errsTrain = errsTest; 
%corsTrain = errsTrain; pcorsTrain = corsTrain;

for s = 1:length(sys)
    for k = 1:length(hzns)
        t1 = hzns(k) + 1;
        %errsTrain(k,:,s) = pRMSE(ypTrain{k,s}(t1:end,:), dataTrain(t1:end,:));
        errsTest(k,:,s) = pRMSE(ypTest{k,s}(t1:end,:), dataTest(t1:end,:));
        for c = 1:width(dataTrain)
            %{
            [corsTrain(k,c,s), pcorsTrain(k,c,s)] = ...
                corr(ypTrain{k,s}{t1:end,c}, dataTrain{t1:end,c});
            %}
            [corsTest(k,c,s), pcorsTest(k,c,s)] = ...
                corr(ypTest{k,s}{t1:end,c}, dataTest{t1:end,c});
        end
    end
end

figStem = figure('Units','normalized', 'Position',[.05,.05,.9,.9]);
K = length(hzns);
for k = 1:K

    subplot(K,3,3*(k-1)+1);
    for s = 1:length(sys)
        colr = sysColr{s};
        stem(errsTest(k,:,s), ['-',colr,'s']); hold on;
    end
    ylabel('pRMSE'); 
    title([num2str(hzns(k)),'-step prediction']);
    xticks(1:width(dataTrain)); xticklabels(dataTrain.Properties.VariableNames);
    xlabel('Channel Name')

    subplot(K,3,3*(k-1)+2);
    for s = 1:length(sys)
        colr = sysColr{s};
        stem(corsTest(k,:,s), ['-',colr,'s']); hold on;
    end
    ylabel('corr'); 
    title([num2str(hzns(k)),'-step prediction']);
    xticks(1:width(dataTrain)); xticklabels(dataTrain.Properties.VariableNames);
    xlabel('Channel Name')

    subplot(K,3,3*(k-1)+3);
    for s = 1:length(sys)
        colr = sysColr{s};
        stem(pcorsTest(k,:,s), ['-',colr,'s']); hold on;
    end
    ylabel('corr p'); 
    title([num2str(hzns(k)),'-step prediction']);
    xticks(1:width(dataTrain)); xticklabels(dataTrain.Properties.VariableNames);
    xlabel('Channel Name')

end

legend(sysName);

%% evaluation - visualization 

% best channel 
%corBest = min(corsTest, [], 3); % best model 
corBest = corsTest(:,:,2); % use bgLTI 
corRepr = abs(corBest - mean(corBest,2));
for k = 1:length(hzns)
    corBest(k,:) = assignrank(corBest(k,:)); % rank channel; higher = better 
    corRepr(k,:) = assignrank(corRepr(k,:)); % lower = more representative 
end
corBest = mean(corBest,1); % mean rank each horizon by channel 
corRepr = mean(corRepr,1);
[~,ch1] = max(corBest); [~,ch2] = min(corRepr);
ch = [ch1; ch2];
H = height(ch);

% plotting
figTime = figure('Units','normalized', 'Position',[.05,.05,.9,.9]);
for c = 1:H

    %{
    % train - timeseries 
    ax(c,1) = subplot(H,4, 4*(c-1)+ 1);
    plottbl(dataTrain, ch(c), ':k', 4); grid on; hold on;
    for s = 1:width(ypTrain)
        for k = 1:length(hzns)
            plottbl(ypTrain{k,s}, ch(c), sysColr{s}, hznlwd(k)); 
        end
    end
    title('Train');

    % train - scatter 
    ax(c,2) = subplot(H,4, 4*(c-1)+ 2); hold on;
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
    ax(c,1) = subplot(H,2, 2*(c-1)+ 1);
    plottbl(dataTest, ch(c), ':k', 4); grid on; hold on;
    for s = 1:width(ypTest)
        for k = 1:length(hzns)
            plottbl(ypTest{k,s}, ch(c), sysColr{s}, hznlwd(k)); 
        end
    end
    title('Test');

    % test - scatter 
    ax(c,2) = subplot(H,2, 2*(c-1)+ 2); hold on;
    for s = 1:width(ypTest)
        colr = sysColr{s};
        for k = 1:length(hzns)
            mkr = hznmkr{k};
            plot(dataTest{:,ch(c)}, ypTrain{k,s}{:,ch(c)}, [mkr,colr]);
        end
    end
    grid on; xlabel('true'); ylabel('pred');
    title('Test');

end

linkaxes(ax(:,1), 'x'); %linkaxes(ax(:,3), 'x');

%% saving 
svname = inputdlg('Save systems as:', 'File Save Name', 1, ...
    {[fn,'_eval']});
if ~isempty(svname)
    svname = svname{1};
    save(fullfile(fp,[svname,'.mat']), 'sys', 'hzns', 'sysName', ...
        'dataTrain', 'dataTest', 'ypTest', 'ypTrain', 'fn', ...
        'corsTest', 'corsTrain', 'pcorsTest', 'pcorsTrain', 'errsTest', 'errsTrain')
    saveas(figStem, fullfile(fp,[svname,'-stem']),'fig'); 
    saveas(figStem, fullfile(fp,[svname,'-stem']),'png'); 
    saveas(figTime, fullfile(fp,[svname,'-time']),'fig'); 
    saveas(figTime, fullfile(fp,[svname,'-time']),'png'); 
end

%% helpers

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