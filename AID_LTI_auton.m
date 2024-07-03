function [trainPred, testPred, trainEval, testEval, A_t_train, A_t_test] = ...
    AID_LTI_auton(trainData, testData, Am, KA, shutoff)

% handle incomplete args 
if nargin < 4
    KA = [];
    if nargin < 3
        Am = [];
        if nargin < 2
            testData = trainData;
            trainData = [];
        end
    end
end
if nargin < 5
    shutoff = false(1,height(testData));
end

% handle empty args 
if isempty(Am)
    Am = -eye(width(testData));
end
if isempty(KA)
    KA = lyap(Am',eye(size(Am)));
end

% check Am Hurwitz 
lambda = eig(Am); 
if sum(lambda >= 0)
    error('Am should be Hurwitz.')
end
% check KA PDS 
lambda = eig(KA);
if sum(lambda <= 0) | ~isequal(KA,KA')
    error('KA should be PDS.')
end
% check lyap deriv NDS
dLyap = Am'*KA + KA*Am; 
lambda = eig(-dLyap);
if sum(lambda <= 0) | ~isequal(dLyap,dLyap')
    error('KA produces unstable system.')
end

%% sample rate 
%{
Time = testData.Time; 
Th = mean(diff(Time));
Th = seconds(Th); % seconds 
%}
Fs = testData.Properties.CustomProperties.SampleRateHz;
Th = 1/Fs; % seconds

%% starting estimate of discrete A 
if isempty(trainData)
    % no starting knowledge 
    A0 = zeros(width(testData));
else
    Xtrain = table2array(trainData)'; 
    X1 = Xtrain(:,1:(end-1)); X2 = Xtrain(:,2:end); 
    A0 = X2*X1' * (X1*X1')^-1;
end

%% run "test" data

Xtest = table2array(testData)';
xtest_est_t = nan(size(Xtest));
%{
A_t_test = nan([size(A0),width(Xtest)]);
A_t_test(:,:,1) = A0;
%}

x = Xtest(:,1); 
A = A0;
xest = x; 
xtest_est_t(:,1) = xest;
[Ad,Bd0] = c2d(Am,eye(size(Am)),Th);
for t = 2:width(Xtest)
    dA = -KA*(xest-x)*x';
    if ~shutoff(t)
        A = A + dA*Th; 
        %[Ad,Bd] = c2d(Am,A-Am,Th);
        Bd = Bd0*(A-Am);
        xest = Ad*xest + Bd*x;
        %dxest = Am*xest + (A-Am)*x;
        %xest = xest + dxest*Th;
    else
        xest = Ad*xest;
    end
    xtest_est_t(:,t) = xest;
    x = Xtest(:,t);
    %A_t_test(:,:,t) = A;
end

testPred = testData; 
testPred{:,:} = xtest_est_t';

%% run "train" data
if ~isempty(trainData)

xtrain_est_t = nan(size(Xtrain));
%{
A_t_train = nan([size(A0),width(Xtrain)]);
A_t_train(:,:,1) = A0;
%}

x = Xtrain(:,1); 
A = A0;
xest = x; 
xtrain_est_t(:,1) = xest;
[Ad,Bd0] = c2d(Am,eye(size(Am)),Th);
for t = 2:width(Xtrain)
    dA = -KA*(xest-x)*x';
        A = A + dA*Th; 
        %[Ad,Bd] = c2d(Am,A-Am,Th);
        Bd = Bd0*(A-Am);
        xest = Ad*xest + Bd*x;
    xtrain_est_t(:,t) = xest;
    x = Xtrain(:,t);
    %A_t_train(:,:,t) = A;
end

trainPred = trainData; 
trainPred{:,:} = xtrain_est_t';

else
    trainPred = [];
end

%% evaluate results 
testEval.RMSE = rmse(Xtest(:), xtest_est_t(:));
testEval.pRMSE = (testEval.RMSE)/rms(Xtest(:));

if ~isempty(trainData)
    trainEval.RMSE = rmse(Xtrain(:), xtrain_est_t(:));
    trainEval.pRMSE = (trainEval.RMSE)/rms(Xtrain(:));
else
    trainEval = [];
end

end