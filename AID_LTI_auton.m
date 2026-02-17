function [trainPred, testPred, trainEval, testEval, A_t_train, A_t_test] = ...
    AID_LTI_auton(trainData, testData, Am, KA, shutoff, showfit)

% handle incomplete args 
if nargin < 6
    showfit = false;
    if nargin < 5
        shutoff = [];
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
end

% handle empty args 
if isempty(Am)
    Am = -eye(width(testData));
elseif numel(Am) == 1
    Am = Am*eye(width(testData));
end
if isempty(KA)
    KA = lyap(Am',eye(size(Am)));
elseif numel(KA) == 1
    KA = lyap(Am',KA*eye(size(Am)));
end
if isempty(shutoff)
    shutoff = false(1,height(testData));
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

AmInv = Am^-1;

% choose display channel 
chandisp = [];
if numel(showfit) > 1
    if iscell(showfit)
        if ischar(showfit{1})
            chandisp = showfit; showfit = true;
        end
    end
end
if isstring(showfit) | ischar(showfit)
    chandisp(showfit); showfit = true;
end
if showfit & isempty(chandisp)
    % choose the most average channel
    xavg = mean(testData.Variables,2);
    xdev = testData.Variables - xavg;
    xdev = sum(xdev.^2);
    [~,ch] = min(xdev);
    chans = testData.Properties.VariableNames;
    if isstring(chans)
        chandisp = chans(ch);
    else
        chandisp = chans{ch};
    end
end

%% sample rate 
%{
Time = testData.Time; 
Th = mean(diff(Time));
Th = seconds(Th); % seconds 
%}
%{
Fs = testData.Properties.SampleRate;
Th = 1/Fs; % seconds
%}
Th = seconds(testData.Properties.TimeStep);

%% starting estimate of discrete A 
if isempty(trainData)
    % no starting knowledge 
    A0 = eye(width(testData)); % discrete 
else
    Xtrain = table2array(trainData); 
    X1 = Xtrain(1:(end-1),:); X2 = Xtrain(2:end,:); 
    Xtrain = Xtrain';
    A0 = (X1\X2)';
end

%% run "test" data

Xtest = table2array(testData)';
xtest_est_t = nan(size(Xtest));

x = Xtest(:,1); 
Ahat_d = A0;
xest = x; 
xtest_est_t(:,1) = xest;
A_d = expm(Am*Th);
Beta_d = Ahat_d - A_d;
M = Th*AmInv*(A_d-eye(size(A_d)));

%%{
A_t_test = nan([size(A0),width(Xtest)]);
A_t_test(:,:,1) = Ahat_d;
%}

if showfit 
    figure('Units','Normalized','Position',[.1 .3 .8 .4]); 
    subplot(121); imgL = imagesc(A0); colorbar; 
    title('Disc A LSQ fit');
    cmin = min(A0(:)); cmax = max(A0(:)); cdiff = cmax-cmin; 
    if cdiff == 0
        cdiff = .01;
    end
    cmin = cmin - .05*cdiff; cmax = cmax + .05*cdiff;
    imgL.Parent.CLim = [cmin, cmax];
    subplot(122); img = imagesc(Ahat_d, [cmin,cmax]); colorbar;
    ttxt = title('Disc A Adapt: 0%');
    pause(1e-9);
    prog_curr = 0; prog_prev = 0; prog_updTick = 1; % percent progress
end

for t = 2:width(Xtest)
    dAhat = -KA*(xest-x)*x';
    delBeta_d = M*dAhat;
    if ~shutoff(t)
        % adaptively update the system 
        Beta_d = Beta_d + delBeta_d;
        xest = A_d*xest + Beta_d*x;
    else
        %%{
        if ~shutoff(t-1)
            % at start of shutoff, should xest be manually set to last reliable x??
            xest = x;
        end
        %}
        % no adaptive updates during shutoff
        xest = A_d*xest + Beta_d*x;
    end
    xtest_est_t(:,t) = xest;
    x = Xtest(:,t);
    Ahat_d = Beta_d + A_d;
    A_t_test(:,:,t) = Ahat_d;
    if showfit
        prog_curr = 100*t/width(Xtest); 
        if prog_curr - prog_prev > prog_updTick
            prog_prev = prog_curr;
            [cmin,chg_min] = min([cmin, min(Ahat_d(:))]); chg_min = chg_min-1;
            [cmax,chg_max] = max([cmax, max(Ahat_d(:))]); chg_max = chg_max-1;
            img.CData = Ahat_d;
            if chg_max | chg_min
                img.Parent.CLim = [cmin, cmax];
                imgL.Parent.CLim = [cmin, cmax];
            end
            ttxt.String = ['Disc A Estimate: ',num2str(prog_curr,3),'%'];
            pause(1e-9);
        end
    end
end

testPred = testData; 
testPred{:,:} = xtest_est_t';

if showfit
    figure; 
    ax(1) = subplot(2,1,1);
    plot(testData, chandisp); 
    grid on; hold on;
    plot(testPred, chandisp);
    legend('Actual', 'Estimate');
    title('Test Data');
    ax(2) = subplot(2,1,2);
    semilogy(testData.Time, (testPred.(chandisp)-testData.(chandisp)).^2);
    grid on;
    title('Squared Error');
    linkaxes(ax, 'x');
end

%% run "train" data
if ~isempty(trainData)

xtrain_est_t = nan(size(Xtrain));
x = Xtrain(:,1); 
xest = x; 
xtrain_est_t(:,1) = xest;
Ahat_d = eye(size(A0));
Beta_d = Ahat_d - A_d;

%%{
A_t_train = nan([size(A0),width(Xtrain)]);
A_t_train(:,:,1) = Ahat_d;
%}

if showfit 
    figure('Units','Normalized','Position',[.1 .3 .8 .4]); 
    subplot(121); imgL = imagesc(A0); colorbar; 
    title('Disc A LSQ fit');
    cmin = min(A0(:)); cmax = max(A0(:)); cdiff = cmax-cmin; 
    cmin = cmin - .05*cdiff; cmax = cmax + .05*cdiff;
    imgL.Parent.CLim = [cmin, cmax];
    subplot(122); img = imagesc(Ahat_d, [cmin,cmax]); colorbar;
    ttxt = title('Disc A Adapt: 0%');
    pause(1e-9);
    prog_curr = 0; prog_prev = 0; prog_updTick = 1; % percent progress
end

for t = 2:width(Xtrain)
    dAhat = -KA*(xest-x)*x';
    delBeta_d = M*dAhat;
        Beta_d = Beta_d + delBeta_d;
        xest = A_d*xest + Beta_d*x;
    xtrain_est_t(:,t) = xest;
    x = Xtrain(:,t);
    Ahat_d = Beta_d + A_d;
    A_t_train(:,:,t) = Ahat_d;

    if showfit
        prog_curr = 100*t/width(Xtrain); 
        if prog_curr - prog_prev > prog_updTick
            prog_prev = prog_curr;
            [cmin,chg_min] = min([cmin, min(Ahat_d(:))]); chg_min = chg_min-1;
            [cmax,chg_max] = max([cmax, max(Ahat_d(:))]); chg_max = chg_max-1;
            img.CData = Ahat_d;
            if chg_max | chg_min
                img.Parent.CLim = [cmin, cmax];
                imgL.Parent.CLim = [cmin, cmax];
            end
            ttxt.String = ['Disc A Estimate: ',num2str(prog_curr,3),'%'];
            pause(1e-9);
        end
    end
end

trainPred = trainData; 
trainPred{:,:} = xtrain_est_t';

if showfit
    figure; 
    ax(1) = subplot(2,1,1);
    plot(trainData, chandisp); 
    grid on; hold on;
    plot(trainPred, chandisp);
    legend('Actual', 'Estimate');
    title('Test Data');
    ax(2) = subplot(2,1,2);
    semilogy(trainData.Time, (trainPred.(chandisp)-trainData.(chandisp)).^2);
    grid on;
    title('Squared Error');
    linkaxes(ax, 'x');
end

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
    trainEval = []; A_t_train = [];
end

end