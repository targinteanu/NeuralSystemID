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
Ac = Ac0;

if showfit 
    figure('Units','Normalized','Position',[.1 .3 .8 .4]); 
    subplot(121); imgL = imagesc(Ac); colorbar; 
    title('Cont A LSQ fit');
    cmin = min(Ac(:)); cmax = max(Ac(:)); cdiff = cmax-cmin; 
    cmin = cmin - .05*cdiff; cmax = cmax + .05*cdiff;
    imgL.Parent.CLim = [cmin, cmax];
    Ac = zeros(size(Ac)); % reset 
    subplot(122); img = imagesc(Ac, [cmin,cmax]); colorbar;
    ttxt = title('Cont A Adapt: 0%');
    pause(1e-9);
    prog_curr = 0; prog_prev = 0; prog_updTick = 1; % percent progress
end

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
        Ad2 = expm(Ac*Th);
        xest = Ad2*xest; % at start of shutoff, should xest be manually set to last reliable x??
    end
    xtest_est_t(:,t) = xest;
    x = Xtest(:,t);
    %A_t_test(:,:,t) = Ac;
    if showfit
        prog_curr = 100*t/width(Xtest); 
        if prog_curr - prog_prev > prog_updTick
            prog_prev = prog_curr;
            [cmin,chg_min] = min([cmin, min(Ac(:))]); chg_min = chg_min-1;
            [cmax,chg_max] = max([cmax, max(Ac(:))]); chg_max = chg_max-1;
            img.CData = Ac;
            if chg_max | chg_min
                img.Parent.CLim = [cmin, cmax];
                imgL.Parent.CLim = [cmin, cmax];
            end
            ttxt.String = ['Cont A Estimate: ',num2str(prog_curr,3),'%'];
            pause(1e-9);
        end
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

if showfit 
    figure('Units','Normalized','Position',[.1 .3 .8 .4]); 
    subplot(121); imgL = imagesc(Ac); colorbar; 
    title('Cont A LSQ fit');
    cmin = min(Ac(:)); cmax = max(Ac(:)); cdiff = cmax-cmin; 
    cmin = cmin - .05*cdiff; cmax = cmax + .05*cdiff;
    imgL.Parent.CLim = [cmin, cmax];
    Ac = zeros(size(Ac)); % reset 
    subplot(122); img = imagesc(Ac, [cmin,cmax]); colorbar;
    ttxt = title('Cont A Adapt: 0%');
    pause(1e-9);
    prog_curr = 0; prog_prev = 0; prog_updTick = 1; % percent progress
end

for t = 2:width(Xtrain)
    dA = -KA*(xest-x)*x';
        Ac = Ac + dA*Th; 
        %[Ad,Bd] = c2d(Am,Ac-Am,Th);
        Bd = Bd0*(Ac-Am);
        xest = Ad*xest + Bd*x;
    xtrain_est_t(:,t) = xest;
    x = Xtrain(:,t);
    %A_t_train(:,:,t) = Ac;
    if showfit
        prog_curr = 100*t/width(Xtrain); 
        if prog_curr - prog_prev > prog_updTick
            prog_prev = prog_curr;
            [cmin,chg_min] = min([cmin, min(Ac(:))]); chg_min = chg_min-1;
            [cmax,chg_max] = max([cmax, max(Ac(:))]); chg_max = chg_max-1;
            img.CData = Ac;
            if chg_max | chg_min
                img.Parent.CLim = [cmin, cmax];
                imgL.Parent.CLim = [cmin, cmax];
            end
            ttxt.String = ['Cont A Estimate: ',num2str(prog_curr,3),'%'];
            pause(1e-9);
        end
    end
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