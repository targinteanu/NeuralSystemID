function [trainPred, testPred, trainEval, testEval, A_t] = ...
    AID_LTI_auton(trainData, testData, Am, KA)

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

% handle empty args 
if isempty(KA)
    KA = eye(width(testData));
end
if isempty(Am)
    Am = -eye(width(testData));
end

% check Am Hurwitz 
lambda = eig(Am); 
if sum(lambda >= 0)
    error('Am should be Hurwitz.')
end
% check KA PDS 
lambda = eig(KA);
if sum(lambda <= 0) | ~equals(KA,KA')
    error('KA should be PDS.')
end

% starting estimate of discrete A 
if isempty(trainData)
    % no starting knowledge 
    A = zeros(width(testData));
else
    Xtrain = trainData{:,:}'; 
    X1 = Xtrain(:,1:(end-1)); X2 = Xtrain(:,2:end); 
    A = X2*X1' * (X1*X1')^-1;
end


end