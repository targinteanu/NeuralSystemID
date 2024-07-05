function [trainPred, testPred, trainEval, testEval, A_t_train, A_t_test] = ...
    AID_LTI_auton(trainData, testData, Am, KA, shutoff, showfit)

% handle incomplete args 
if nargin < 5
    showfit = false;
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
    warning('Am should be Hurwitz.')
end
% check KA PDS 
lambda = eig(KA);
if sum(lambda <= 0) | ~isequal(KA,KA')
    warning('KA should be PDS.')
end
% check lyap deriv NDS
dLyap = Am'*KA + KA*Am; 
lambda = eig(-dLyap);
if sum(lambda <= 0) | ~isequal(dLyap,dLyap')
    warning('KA produces unstable system.')
end

%% sample rate 
%{
Time = testData.Time; 
Th = mean(diff(Time));
Th = seconds(Th); % seconds 
%}
Fs = testData.Properties.SampleRate;
Th = 1/Fs; % seconds

%% starting estimate of discrete A 
if isempty(trainData)
    % no starting knowledge 
    A0 = eye(width(testData)); % discrete 
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
A_t_test(:,:,1) = Ac0;
%}

x = Xtest(:,1); 
A = A0;
xest = x; 
xtest_est_t(:,1) = xest;
[Ad,Bd0] = c2d(Am,eye(size(Am)),Th);
Ac0 = d2c(A,zeros(size(A,1),1),Th);
if showfit 
    figure('Units','Normalized','Position',[.1 .3 .8 .4]); 
    subplot(121); imagesc(Ac0); colorbar; 
    title('Cont A LSQ fit');
    Ac0 = zeros(size(Ac0)); % reset 
    subplot(122); img = imagesc(Ac0); colorbar;
    ttxt = title('Cont A Adapt: 0%');
    pause(1e-9);
end
Ac = Ac0;
for t = 2:width(Xtest)
    dA = -KA*(xest-x)*x';
    if ~shutoff(t)
        Ac = Ac + dA*Th; 
        %[Ad,Bd] = c2d(Am,Ac-Am,Th);
        Bd = Bd0*(Ac-Am);
        xest = Ad*xest + Bd*x;
        %dxest = Am*xest + (Ac-Am)*x;
        %xest = xest + dxest*Th;
    else
        xest = Bd*xest; % at start of shutoff, should xest be manually set to last reliable x??
    end
    xtest_est_t(:,t) = xest;
    x = Xtest(:,t);
    %A_t_test(:,:,t) = Ac;
    if showfit
        img.CData = Ac;
        ttxt.String = ['Cont A Estimate: ',num2str(100*t/width(Xtest),3),'%'];
        pause(1e-9);
    end
end

testPred = testData; 
testPred{:,:} = xtest_est_t';

%% run "train" data
if ~isempty(trainData)

xtrain_est_t = nan(size(Xtrain));
%{
A_t_train = nan([size(A0),width(Xtrain)]);
A_t_train(:,:,1) = Ac0;
%}

x = Xtrain(:,1); 
A = A0;
xest = x; 
xtrain_est_t(:,1) = xest;
Ac = Ac0;
for t = 2:width(Xtrain)
    dA = -KA*(xest-x)*x';
        Ac = Ac + dA*Th; 
        %[Ad,Bd] = c2d(Am,Ac-Am,Th);
        Bd = Bd0*(Ac-Am);
        xest = Ad*xest + Bd*x;
    xtrain_est_t(:,t) = xest;
    x = Xtrain(:,t);
    %A_t_train(:,:,t) = Ac;
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