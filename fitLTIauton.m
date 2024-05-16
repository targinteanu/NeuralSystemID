function [trainPred, testPred] = fitLTIauton(xData)

X = xData{:,:}'; 
X1 = X(:,1:(end-1)); X2 = X(:,2:end); 
A = X2*X1' * (X1*X1')^-1;

trainPred = xData; 
for t = 2:height(trainPred)
    %x0 = trainPred{t-1,:}'; 
    x0 = xData{t-1,:}';
    x = A*x0; 
    trainPred{t,:} = x';
end

testPred = xData;
%testPred{1,:} = x';
testPred{1,:} = xData{end,:};
for t = 2:height(testPred)
    x0 = testPred{t-1,:}'; 
    x = A*x0; 
    testPred{t,:} = x';
end
t = testPred.Time; 
t = t - t(1) + trainPred.Time(end); 
testPred.Time = t;

end