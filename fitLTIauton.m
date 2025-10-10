function [trainPred, testPred, trainEval, testEval, A] = ...
    fitLTIauton(trainData, testData, PredStateVisible)

% handle input args 
if nargin < 2
    testData = [];
end
if nargin < 3
    PredStateVisible = true;
end

% NEED TO ACCOUNT FOR DISCONTINUOUS DATA!

% Train A matrix
Xtrain = table2array(trainData)'; 
X1 = Xtrain(:,1:(end-1)); X2 = Xtrain(:,2:end); 
A = X2*X1' * (X1*X1')^-1;

% setup prediction outputs 
trainPredX = Xtrain; 
for t = 2:width(trainPredX)
    if PredStateVisible
        x0 = Xtrain(:,t-1);
    else
        x0 = trainPredX(:,t-1); 
    end
    x = A*x0; 
    trainPredX(:,t) = x;
end

% compare prediction to testing data 
if ~isempty(testData) % consider consolidating by replacing this with projLTIauton
Xtest = table2array(testData)';
testPredX = Xtest;
%testPredX = Xtrain;
%testPredX(:,1) = x;
%testPredX(:,1) = Xtrain(:,end);
for t = 2:width(testPredX)
    if PredStateVisible
        x0 = Xtest(:,t-1);
    else
        x0 = testPredX(:,t-1); 
    end
    x = A*x0; 
    testPredX(:,t) = x;
end
testPred = testData; 
testPred{:,:} = testPredX';
%t = testPred.Time; 
%t = t - t(1) + trainPred.Time(end); 
%testPred.Time = t;
end

% format outputs 
trainPred = trainData; 
trainPred{:,:} = trainPredX';

% output evaluation stats - training accuracy 
trainEval.RMSE = rmse(Xtrain(:), trainPredX(:)); 
trainEval.pRMSE = rmse(Xtrain(:), trainPredX(:))/rms(Xtrain(:));
%{
[r,p] = corr(Xtrain', trainPredX'); 
r = diag(r); p = diag(p); 
trainEval.CorrMean = mean(r); trainEval.CorrStd = std(r); 
trainEval.CorrP = 1-prod(1-p);
%}

% output evaluation stats - testing accuracy
if ~isempty(testData)
testEval.RMSE = rmse(Xtest(:), testPredX(:)); 
testEval.pRMSE = rmse(Xtest(:), testPredX(:))/rms(Xtest(:));
%{
[r,p] = corr(Xtest', testPredX'); 
r = diag(r); p = diag(p); 
testEval.CorrMean = mean(r); testEval.CorrStd = std(r); 
testEval.CorrP = 1-prod(1-p);
%}
else
    testPred = [];
    testEval.RMSE = nan;
    testEval.pRMSE = nan;
end

end