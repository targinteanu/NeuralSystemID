%% Parkinson's Disease (PD) Project - evaluate autonomous systems
% Autonomous only: does not include brain stimulation. 

%% load data file 
[fn,fp] = uigetfile('trainednet*andsys*.mat');
load(fullfile(fp,fn), 'bgNSS', 'sysNull', 'sysLTI', 'sysAR', 'bgLTI', ...
    'dataTrain', 'dataTest');
disp([fp,' --- ',fn]);
[~,fn] = fileparts(fn);
sysName = {'bgNSS', 'sysNull', 'sysLTI', 'sysAR', 'bgLTI'};
sysColr = {'-b',    '-k',      '-r',     '-c',    '-m'};
sys = cellfun(@eval,sysName, 'UniformOutput',false);
fsNew = dataTrain.Properties.SampleRate;

%% select training/testing validation (subsets) 
Lval = 5000; % # samples 
disp([num2str(Lval/fsNew),' seconds selected for validation.'])
[~,dataTestSel] = max(cellfun(@height, dataTest));
dataTest = dataTest{dataTestSel};
dataTest = dataTest(1:Lval,:); dataTrain = dataTrain(1:Lval,:);

%% k-step ahead prediction 
% artifact duration is around 10ms, but we should verify different time
% scales. 
hzns = [.025, .25, 2.5]; % seconds 
hzns = ceil(hzns * fsNew); % # samples 

ypTrain = {}; ypTest = {};
for hzn = hzns
    disp(['Eval ',num2str(hzn),'-step']);
    disp('Training Eval');
    ypTrain_ = cellfun(@(S) myPredict(S, dataTrain, hzn, true), sys, ...
        'UniformOutput',false);
    ypTrain = [ypTrain; ypTrain_];
    disp('Testing Eval');
    ypTest_ = cellfun(@(S) myPredict(S, dataTest, hzn, true), sys, ...
        'UniformOutput',false);
    ypTest = [ypTest; ypTest_];
end
clear ypTrain_ ypTest_ 

%% evaluation - assessment

pRMSE = @(yp, y) rmse(yp.Variables, y.Variables) ./ rms(y.Variables);
rho = @(yp, y) arrayfun(@(c) corr(y{:,c}, yp{:,c}), 1:width(y));

errsTrain = nan(length(hzns), width(dataTrain), length(sys)); 
corsTrain = errsTrain;
errsTest = errsTrain; 
corsTest = errsTest; 

for s = 1:length(sys)
    for k = 1:length(hzns)
        errsTrain(k,:,s) = pRMSE(ypTrain{k,s}, dataTrain);
        corsTrain(k,:,s) = rho(ypTrain{k,s}, dataTrain);
        errsTest(k,:,s) = pRMSE(ypTest{k,s}, dataTest);
        corsTest(k,:,s) = rho(ypTest{k,s}, dataTest);
    end
end

figStem = figure('Units','normalized', 'Position',[.05,.05,.9,.9]);
K = length(hzns);
for k = 1:K

    subplot(K,2,2*(k-1)+1);
    for s = 1:length(sys)
        colr = sysColr{s};
        stem(errsTrain(k,:,s), [colr,'s']); hold on;
        colr(1) = ':';
        stem(errsTest(k,:,s), [colr,'x'], 'LineWidth',1);
    end
    ylabel('pRMSE'); 
    title([num2str(hzns(k)),'-step prediction']);
    xticks(1:width(dataTrain)); xticklabels(dataTrain.Properties.VariableNames);
    xlabel('Channel Name')

    subplot(K,2,2*(k-1)+2);
    for s = 1:length(sys)
        colr = sysColr{s};
        stem(corsTrain(k,:,s), [colr,'s']); hold on;
        colr(1) = ':';
        stem(corsTest(k,:,s), [colr,'x'], 'LineWidth',1);
    end
    ylabel('corr'); 
    title([num2str(hzns(k)),'-step prediction']);
    xticks(1:width(dataTrain)); xticklabels(dataTrain.Properties.VariableNames);
    xlabel('Channel Name')

end

lgd = [sysName; sysName];
legend(lgd(:));

%% evaluation - visualization 

% best channel 
%corBest = min(corsTest, [], 3); % best model 
corBest = corsTest(:,:,5); % use bgLTI 
corRepr = abs(corBest - mean(corBest,2));
for k = 1:length(hzns)
    corBest(k,:) = assignrank(corBest(k,:)); % rank channel; higher = better 
    corRepr(k,:) = assignrank(corRepr(k,:)); % lower = more representative 
end
corBest = mean(corBest,1); % mean rank each horizon by channel 
corRepr = mean(corRepr,1);
[~,ch1] = max(corBest); [~,ch2] = min(corRepr);

% get at least one of each preprocessing type (power, freq)
ch3 = mod(ch1+width(dataTrain)/2, width(dataTrain));
ch4 = mod(ch2+width(dataTrain)/2, width(dataTrain));
ch = [ch1; ch3; ch2; ch4]; ch = round(ch);
H = height(ch);

% plotting
figTime = figure('Units','normalized', 'Position',[.05,.05,.9,.9]);
for c = 1:H

    % train - timeseries 
    ax(c,1) = subplot(H,4, 4*(c-1)+ 1);
    plottbl(dataTrain, ch(c), ':k', 4); grid on; hold on;
    for s = 1:width(ypTrain)
        plottbl(ypTrain{1,s}, ch(c), sysColr{s}, 2); 
        plottbl(ypTrain{2,s}, ch(c), sysColr{s}, 1.25); 
        plottbl(ypTrain{3,s}, ch(c), sysColr{s}, .5); 
    end
    title('Train');

    % train - scatter 
    ax(c,2) = subplot(H,4, 4*(c-1)+ 2); hold on;
    for s = 1:width(ypTrain)
        colr = sysColr{s}(2);
        plot(dataTrain{:,ch(c)}, ypTrain{1,s}{:,ch(c)}, ['d',colr]);
        plot(dataTrain{:,ch(c)}, ypTrain{2,s}{:,ch(c)}, ['o',colr]);
        plot(dataTrain{:,ch(c)}, ypTrain{3,s}{:,ch(c)}, ['.',colr]);
    end
    grid on; xlabel('true'); ylabel('pred');
    title('Train');

    % test - timeseries 
    ax(c,3) = subplot(H,4, 4*(c-1)+ 3);
    plottbl(dataTest, ch(c), ':k', 4); grid on; hold on;
    for s = 1:width(ypTest)
        plottbl(ypTest{1,s}, ch(c), sysColr{s}, 2); 
        plottbl(ypTest{2,s}, ch(c), sysColr{s}, 1.25); 
        plottbl(ypTest{3,s}, ch(c), sysColr{s}, .5); 
    end
    title('Test');

    % test - scatter 
    ax(c,4) = subplot(H,4, 4*(c-1)+ 4); hold on;
    for s = 1:width(ypTest)
        colr = sysColr{s}(2);
        plot(dataTest{:,ch(c)}, ypTest{1,s}{:,ch(c)}, ['d',colr]);
        plot(dataTest{:,ch(c)}, ypTest{2,s}{:,ch(c)}, ['o',colr]);
        plot(dataTest{:,ch(c)}, ypTest{3,s}{:,ch(c)}, ['.',colr]);
    end
    grid on; xlabel('true'); ylabel('pred');
    title('Test');

end

linkaxes(ax(:,1), 'x'); linkaxes(ax(:,3), 'x');

%% saving 
svname = inputdlg('Save systems as:', 'File Save Name', 1, ...
    {[fn,'_eval']});
if ~isempty(svname)
    svname = svname{1};
    save(fullfile(fp,[svname,'.mat']), 'sys', 'hzns', 'sysColr', 'sysName', ...
        'dataTrain', 'dataTest', 'ypTest', 'ypTrain', 'fn')
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