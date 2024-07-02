function [trainPred, testPred, trainEval, testEval, A] = ...
    fitLTIauton(trainData, testData)

if nargin < 2
    testData = [];
end

Xtrain = table2array(trainData)'; 
X1 = Xtrain(:,1:(end-1)); X2 = Xtrain(:,2:end); 
A = X2*X1' * (X1*X1')^-1;

trainPred = trainData; 
for t = 2:height(trainPred)
    %x0 = table2array(trainPred(t-1,:))'; 
    x0 = table2array(trainData(t-1,:))';
    x = A*x0; 
    trainPred{t,:} = x';
end

if ~isempty(testData)
%testPred = trainData;
testPred = testData;
%testPred{1,:} = x';
%testPred{1,:} = table2array(trainData(end,:));
for t = 2:height(testPred)
    %x0 = table2array(testPred(t-1,:))'; 
    x0 = table2array(testData(t-1,:))';
    x = A*x0; 
    testPred{t,:} = x';
end
%t = testPred.Time; 
%t = t - t(1) + trainPred.Time(end); 
%testPred.Time = t;
testPredX = table2array(testPred(:,:))'; 
Xtest = table2array(testData(:,:))'; 
end

trainPredX = table2array(trainPred(:,:))'; 

trainEval.RMSE = rmse(Xtrain(:), trainPredX(:)); 
trainEval.pRMSE = rmse(Xtrain(:), trainPredX(:))/rms(Xtrain(:));
%{
[r,p] = corr(Xtrain', trainPredX'); 
r = diag(r); p = diag(p); 
trainEval.CorrMean = mean(r); trainEval.CorrStd = std(r); 
trainEval.CorrP = 1-prod(1-p);
%}

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