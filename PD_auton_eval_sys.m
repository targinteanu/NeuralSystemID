%% Parkinson's Disease (PD) Project - evaluate autonomous systems
% Autonomous only: does not include brain stimulation. 

%% load data file 
[fn,fp] = uigetfile('trainednet*andsys*.mat');
load(fullfile(fp,fn), 'bgNSS', 'sysNull', 'sysLTI', 'sysAR', 'bgLTI', ...
    'dataTrain', 'dataTest');
disp([fp,' --- ',fn]);
sysName = {'bgNSS', 'sysNull', 'sysLTI', 'sysAR', 'bgLTI'};
sysColr = {'-b',    '-k',      '-r',     '-c',    '-m'};
sys = cellfun(@eval,sysName, 'UniformOutput',false);

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
    ypTest_ = {};
    for trl = 1:height(dataTest)
        disp(['Testing Eval ',num2str(trl)]);
        dataTest_ = dataTest{trl};
        ypTest__ = cellfun(@(S) myPredict(S, dataTest_, hzn, true), sys, ...
            'UniformOutput',false);
        ypTest_ = [ypTest_; ypTest__];
        if height(ypTest_) > 1
            for s = 1:width(ypTest_)
                ypTest_{1,s} = [ypTest_{1,s}; ypTest_{2,s}];
            end
            ypTest_ = ypTest(1,:);
        end
    end
    ypTest = [ypTest; ypTest_];
end
clear ypTrain_ dataTest_ ypTest_ ypTest__

for trl = 1:height(dataTest)
    dataTest{1,1} = [dataTest{1,1}; dataTest{trl,1}];
end
dataTest = dataTest{1,1};

%% evaluation 

pRMSE = @(yp, y) rmse(yp.Variables, y.Variables) ./ rms(y.Variables);
rho = @(yp, y) arrayfun(@(c) corr(y{:,c}, yp{:,c}), 1:width(y));

errsTrain = nan(length(sys), width(dataTrain), length(hzns)); 
corsTrain = errsTrain;
errsTest = errsTrain; 
corsTest = errsTest; 

for s = 1:length(sys)
    for k = 1:length(hzns)
        errsTrain(s,:,k) = pRMSE(ypTrain{s,k}, dataTrain);
        corsTrain(s,:,k) = rho(ypTrain{s,k}, dataTrain);
        errsTest(s,:,k) = pRMSE(ypTest{s,k}, dataTest);
        corsTest(s,:,k) = rho(ypTest{s,k}, dataTest);
    end
end

figStem = figure('Units','normalized', 'Position',[.1,.1,.8,.8]);
K = length(hzns);
for k = 1:K

    subplot(K,2,2*(k-1)+1);
    for s = 1:length(sys)
        colr = sysColr{s};
        stem(errsTrain(s,:,k), [colr,'s']);
        colr(1) = ':';
        stem(errsTest(s,:,k), [colr,'x'], 'LineWidth',1);
    end
    ylabel('pRMSE'); 
    title([num2str(hzns(k)),'-step prediction']);
    xticks(1:width(dataTrain)); xticklabels(dataTrain.Properties.VariableNames);
    xlabel('Channel Name')

    subplot(K,2,2*(k-1)+2);
    for s = 1:length(sys)
        colr = sysColr{s};
        stem(corsTrain(s,:,k), [colr,'s']);
        colr(1) = ':';
        stem(corsTest(s,:,k), [colr,'x'], 'LineWidth',1);
    end
    ylabel('corr'); 
    title([num2str(hzns(k)),'-step prediction']);
    xticks(1:width(dataTrain)); xticklabels(dataTrain.Properties.VariableNames);
    xlabel('Channel Name')

end