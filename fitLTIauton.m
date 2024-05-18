function [trainPred, testPred, trainEval, testEval, A] = ...
    fitLTIauton(trainData, testData)

Xtrain = trainData{:,:}'; 
X1 = Xtrain(:,1:(end-1)); X2 = Xtrain(:,2:end); 
A = X2*X1' * (X1*X1')^-1;

trainPred = trainData; 
for t = 2:height(trainPred)
    %x0 = trainPred{t-1,:}'; 
    x0 = trainData{t-1,:}';
    x = A*x0; 
    trainPred{t,:} = x';
end

%testPred = trainData;
testPred = testData;
%testPred{1,:} = x';
%testPred{1,:} = trainData{end,:};
for t = 2:height(testPred)
    %x0 = testPred{t-1,:}'; 
    x0 = testData{t-1,:}';
    x = A*x0; 
    testPred{t,:} = x';
end
%t = testPred.Time; 
%t = t - t(1) + trainPred.Time(end); 
%testPred.Time = t;

trainPredX = trainPred{:,:}'; testPredX = testPred{:,:}'; 
Xtest = testData{:,:}'; 

trainEval.RMSE = rmse(Xtrain(:), trainPredX(:)); 
[r,p] = corr(Xtrain', trainPredX'); 
r = diag(r); p = diag(p); 
trainEval.CorrMean = mean(r); trainEval.CorrStd = std(r); 
trainEval.CorrP = 1-prod(1-p);

testEval.RMSE = rmse(Xtest(:), testPredX(:)); 
[r,p] = corr(Xtest', testPredX'); 
r = diag(r); p = diag(p); 
testEval.CorrMean = mean(r); testEval.CorrStd = std(r); 
testEval.CorrP = 1-prod(1-p);

end