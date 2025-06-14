function [testPred, trainPred, trainEval, testEval, W_t_train, W_t_test, ARmdl] = ...
    AdaptAR(trainData, testData, N, stepsize, shutoff, donorm, showfit)

% handle incomplete args 
if nargin < 7
    showfit = false;
    if nargin < 6
        donorm = false;
        if nargin < 5
            shutoff = [];
            if nargin < 4
                stepsize = [];
                if nargin < 3
                    N = [];
                    if nargin < 2
                        testData = trainData;
                        trainData = [];
                    end
                end
            end
        end
    end
end

% handle empty args
if isempty(shutoff)
    shutoff = false(1,height(testData));
end
if isempty(stepsize)
    stepsize = 1e-8;
end
if isempty(N)
    N = 10;
end

% handle multi-variable 
if width(testData) > 1
    error('Testing data should be single variable.')
end
if width(trainData) > 1
    error('Training data should be single variable.')
end

%% starting estimate 
if isempty(trainData)
    % no starting knowledge 
    ARmdl = [];
    W0 = zeros(1,N); 
else
    ARmdl = ar(trainData, N, 'yw'); 
    W0 = ARmdl.A; W01 = W0(1); W0 = W0(2:end); W0 = -W0/W01; % AR model wts
end
W0 = fliplr(W0);

%% run "test" data

W1 = zeros(size(W0)); W_t_test = zeros(height(testData), length(W1));
E = nan(size(testData));
testPred = testData; testPred{:,:} = nan;

if showfit
    figure('Units','Normalized','Position',[.1 .3 .8 .4]); 
    subplot(121); 
    stem(W0); hold on; w1plt = stem(W1); grid on;
    xlabel('tap'); ylabel('weight');
    legend('Yule-Walker', 'Adaptive');
    ttl = title(['Step 0 of ',num2str(height(testData))]);
    pause(.001);
    subplot(122);
    ePlt = plot(E.^2); grid on;
    xlabel('timestep'); ylabel('e^2');
    progticksteps = 1000;
    sgtitle('Testing');
end

% dynamic updating AR by online least squares gradient descent
for t = (N+1):height(testData)
    y = testData{t,1}; 
    x = testData{(t-N):(t-1),1};
    ypred = W1*x; testPred{t,1} = ypred;
    E(t) = y-ypred;
    if ~shutoff(t)
        % update weights only when there is not artifact
        del = x*E(t);
        if donorm
            del = del./(x'*x + eps);
        end
        W1 = W1 + stepsize*del';
    end
    W_t_test(t,:) = W1;
    if showfit
        if ~mod(t,progticksteps)
            w1plt.YData = W1;
            ePlt.YData = E.^2;
            ttl.String = ['Step ',num2str(t),' of ',num2str(height(testData))];
            pause(.001);
        end
    end
end

%% run "train" data
if ~isempty(trainData)

W1 = zeros(size(W0)); W_t_train = zeros(height(trainData), length(W1));
E = nan(size(trainData));
trainPred = trainData; trainPred.Variables = nan;

if showfit
    figure('Units','Normalized','Position',[.1 .3 .8 .4]); 
    subplot(121); 
    stem(W0); hold on; w1plt = stem(W1); grid on;
    xlabel('tap'); ylabel('weight');
    legend('Yule-Walker', 'Adaptive');
    ttl = title(['Step 0 of ',num2str(height(trainData))]);
    pause(.001);
    subplot(122);
    ePlt = plot(E.^2); grid on;
    xlabel('timestep'); ylabel('e^2');
    progticksteps = 1000;
    sgtitle('Training');
end

% dynamic updating AR by online least squares gradient descent
for t = (N+1):height(trainData)
    y = trainData{t,1}; 
    x = trainData{(t-N):(t-1),1};
    ypred = W1*x; trainPred{t,1} = ypred;
    E(t) = y-ypred;
    del = x*E(t);
    if donorm
        del = del./(x'*x + eps);
    end
    W1 = W1 + stepsize*del'; 
    W_t_train(t,:) = W1;
    if showfit
        if ~mod(t,progticksteps)
            w1plt.YData = W1;
            ePlt.YData = E.^2;
            ttl.String = ['Step ',num2str(t),' of ',num2str(height(trainData))];
            pause(.001);
        end
    end
end

else
    trainPred = [];
end

%% evaluate results 
testEval.RMSE = rmse(testData.Variables, testPred.Variables);
testEval.pRMSE = (testEval.RMSE)/rms(testData.Variables);

if ~isempty(trainData)
    trainEval.RMSE = rmse(trainData.Variables, trainPred.Variables);
    trainEval.pRMSE = (trainEval.RMSE)/rms(trainData.Variables);
else
    trainEval = [];
end

end