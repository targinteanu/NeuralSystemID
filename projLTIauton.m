function [testPred,testEval] = projLTIauton(A,testData,shutoff)

if nargin < 3
    shutoff = false(1, height(testData));
end

Xtest = table2array(testData)';
testPredX = Xtest;

for t = 2:width(testPredX)
    if ~shutoff(t)
        x0 = Xtest(:,t-1);
    else
        x0 = testPredX(:,t-1); 
    end
    x = A*x0; 
    testPredX(:,t) = x;
end

testPred = testData; 
testPred{:,:} = testPredX';

testEval.RMSE = rmse(Xtest(:), testPredX(:)); 
testEval.pRMSE = rmse(Xtest(:), testPredX(:))/rms(Xtest(:));

end