function [testPred,testEval] = projAR(mdl,testData,shutoff)

if nargin < 3
    shutoff = false(1,height(testData));
end
if width(testData) > 1
    error('Testing data should be single variable.')
end
if length(shutoff) ~= height(testData)
    error('incompatible input sizes.')
end

% extract inputs 
Ytest = testData.Variables; 
YtestPred = Ytest;
N = length(mdl.A)-1;

% process shutoff start/end
shutoff = [false, shutoff]; % plug in case first timepoint is shut off
shutoffStart = diff(shutoff) > 0; % rising edge 
shutoffEnd = diff(shutoff) < 0; % falling edge
shutoffStart = find(shutoffStart); shutoffEnd = find(shutoffEnd); 
if shutoffStart(end) > shutoffEnd(end)
    % plug in case last timepoint was shut off
    shutoffEnd = [shutoffEnd, length(shutoff)]; 
end
% now shutoffStart and shutoffEnd should be same size 
shutoffStartEnd = [shutoffStart; shutoffEnd];

% main loop 
for idx = 1:length(shutoffStart)

    % "good" forecast before shutoff
    if idx > 1
        t1 = shutoffEnd(idx-1) + 1;
    else
        t1 = 1;
    end
    y = Ytest(1:shutoffStart(idx),:);
    ypred = y;
    for t = (t1+1):height(y)
        yPst = y(1:(t-1),:);
        ypred(t) = myFastForecastWrap(mdl,yPst,1);
    end
    YtestPred((t1+1):shutoffStart(idx),:) = ypred((t1+1):shutoffStart(idx),:);

    % "bad" forecast during shutoff
    yPst = Ytest(1:(shutoffStart(idx)-1),:);
    K = shutoffEnd(idx) - shutoffStart(idx) + 1;
    ypred = myFastForecastWrap(mdl,yPst,K+1);
    YtestPred(shutoffStart(idx):(shutoffEnd(idx)+1),:) = ypred;
    Ytest(shutoffStart(idx):shutoffEnd(idx),:) = ypred(1:K,:);

end

% final "good" forecast
t1 = shutoffEnd(end) + 1;
y = Ytest;
ypred = y;
for t = (t1+1):height(y)
    yPst = y(1:(t-1),:);
    ypred(t) = myFastForecastWrap(mdl,yPst,1);
end
YtestPred((t1+1):end,:) = ypred((t1+1):end,:);

% output assignment 
testPred = testData; 
testPred{:,:} = YtestPred;
testEval.RMSE = rmse(Ytest, YtestPred); 
testEval.pRMSE = rmse(Ytest, YtestPred)/rms(Ytest);

    function yf = myFastForecastWrap(arMdl, yPast, k)
        h = height(yPast);
        if h < N
            yPast = [zeros(N-h,1); yPast];
        end
        yf = myFastForecastAR(arMdl, yPast, k);
    end

end