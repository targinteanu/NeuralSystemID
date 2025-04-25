%% Parkinson's Disease (PD) Project - simulate autonomous systems
% Simulate at every time point in order to figure out how many samples
% ahead the prediction is valid. 
% Autonomous only: does not include brain stimulation. 
% not testing neural network

%% load data file 
[fn,fp] = uigetfile('*andsys*.mat');
load(fullfile(fp,fn), 'sysNull', 'sysLTI', 'sysAR', 'bgLTI', ...
    'dataTrain', 'dataTest');
disp([fp,' --- ',fn]);
[~,fn] = fileparts(fn);
sysName = {'sysAR', 'bgLTI'};
sysColr = {'b',     'r'};
sys = cellfun(@eval,sysName, 'UniformOutput',false);
sysName = {'AR', 'LTI'};
fs = dataTrain.Properties.SampleRate;

%% simulation 
% Ysim will have entries for each system (AR and LTI) and each testing data
% set (should be 2). Entries have timepoints as rows and channels as
% columns. Simulations starting at regular intervals are stacked in the
% third dimension. 

N = width(sysAR.A{1,1}) - 1; % order of AR model
Cinv = pinv(bgLTI.C); % quick convert outputs to estimate of hidden states 
simSkipNum = 10; % new sim starts at every _th sample (for time/memory)
simMaxDur = 5000; % samples; 75 sec
Ysim = cell(height(dataTest), width(sys)); % all sim results 

% prepare starting point data for each set of test data 
for d = 1:height(dataTest)
    dataTest_ = dataTest{d};
    if ~isempty(dataTest_)
        dataTest_ = dataTest_(1:floor(height(dataTest_)/2), :);
        if ~isempty(dataTest_)
            dataTest_ = retime(dataTest_, 'regular','nearest',...
                'TimeStep', seconds(bgLTI.Ts) );
            dataTest_.Time = dataTest_.Time - dataTest_.Time(1);
            t = dataTest_(:,[]); % timetable of time only 
            t = t(1:min(simMaxDur, height(t)), :);

            % LTI simulations 
            sN = ceil(height(dataTest_)/simSkipNum);
            YsLTI = nan(height(t), width(dataTest_), sN );
            si = 1;
            for ti = 1:simSkipNum:height(dataTest_)
                if ~mod(si, 10)
                    disp(['LTI: simulation ',num2str(si),' of ',num2str(sN)])
                end
                x0 = Cinv * dataTest_{ti,:}';
                ysLTI = sim(bgLTI, t, simOptions('InitialCondition',x0));
                YsLTI(:,:,si) = ysLTI.Variables;
                %t.Time = t.Time + t.Properties.TimeStep; % shift to next timestep
                si = si+1;
                clear ysLTI
            end

            % AR simulations 
            sN = ceil( (height(dataTest_)-N)/simSkipNum );
            YsAR = nan(height(t)-N+1, width(dataTest_), sN );
            si = 1;
            for ti = (N+1):simSkipNum:height(dataTest_)
                if ~mod(si, 10)
                    disp(['AR: simulation ',num2str(si),' of ',num2str(sN)])
                end
                ysAR = myFastForecastAR(sysAR, dataTest_{(ti-N):(ti-1),:}, ...
                    height(t)-N);
                ysAR = [dataTest_{ti-1,:}; ysAR]; % ??
                YsAR(:,:,si) = ysAR;
                si = si+1;
                clear ysAR
            end

            Ysim{d,2} = YsLTI; Ysim{d,1} = YsAR;
            %clear YsLTI YsAR
        end
    end
    clear dataTest_
end

%% sim -> pred 
% Ypred entries will have rows corresponding to prediction horizons, and
% these entries should be similar to a downsampled output of MATLAB's
% predict function ("yp") at that horizon, stacked along the third
% dimension with the corresponding actual data ("y").

Ypred = cell(height(dataTest), width(sys)); % test set * sys * hzn

for d = 1:height(dataTest)
    dataTest_ = dataTest{d};
    if ~isempty(dataTest_)
        YsLTI = Ysim{d,2}; YsAR = Ysim{d,1};

        YpLTI = cell(height(YsLTI),1);
        for h = 1:height(YsLTI)
            y = dataTest_{h:simSkipNum:end, :}; 
            %y = y(1:ceil(height(YsLTI)/(simSkipNum)), :); 
            yp = squeeze(YsLTI(h,:,:))';
            y = y(1:height(yp), :);
            YpLTI{h} = single(cat(3,y,yp));
            clear y yp
        end

        YpAR = cell(height(YsAR),1);
        for h = 1:height(YsAR)
            y = dataTest_{(h+N-1):simSkipNum:end, :}; 
            %y = y(1:ceil(height(YsAR)/(simSkipNum)), :); 
            %y = y(N:end,:);
            yp = squeeze(YsAR(h,:,:))';
            y = y(1:height(yp), :);
            YpAR{h} = single(cat(3,y,yp));
            clear y yp
        end

        Ypred{d,2} = YpLTI; Ypred{d,1} = YpAR;
        % clear YpLTI YpAR
    end
end

% collapse test sets 
Ypred_ = cell(1, width(Ypred));
for s = 1:width(Ypred)
    h = cellfun(@height, Ypred(:,s));
    Yp = cell(max(h),1);
    for d = 1:height(Ypred)
        yp = Ypred{d,s};
        for h = 1:height(yp)
            Yp{h} = [Yp{h}; yp{h}];
        end
        clear yp
    end
    Ypred_{s} = Yp;
    clear Yp
end
Ypred = Ypred_; clear Ypred_

%% evaluate pred correlation and error
% In each channel (and across all channels together), evaluate the
% correlation and the fraction RMS error between the predicted and actual
% data as a function of prediction horizon. 

cors = cell(size(Ypred)); pcor = cors; % correlation and p-value 
errs = cell(size(Ypred)); % pRMSE 

for s = 1:length(Ypred)
    disp(['Evaluating ',sysName{s}])
    Yp = Ypred{s};
    cors_ = nan(height(Yp), width(dataTrain)+1); % hzn * chan (incl overall)
    pcor_ = cors_; errs_ = cors_;
    for h = 1:height(cors_) % hzn
        if ~mod(h, 10)
            disp([' - timepoint ',num2str(h),' of ',num2str(height(cors_))])
        end

        % eval each channel
        for c = 1:width(dataTrain)
            y = Yp{h}(:,c,1); yp = Yp{h}(:,c,2);
            [cors_(h,c), pcor_(h,c)] = corr(y,yp,"tail","right");
            errs_(h,c) = rmse(yp,y) ./ rms(y);
        end
        
        % eval overall/agregate channel at end 
        y = Yp{h}(:,:,1); yp = Yp{h}(:,:,2);
        y = y(:); yp = yp(:);
        [cors_(h,c+1), pcor_(h,c+1)] = corr(y,yp,"tail","right");
        errs_(h,c+1) = rmse(yp,y) ./ rms(y);

        clear y yp 
    end
    cors{s} = cors_; pcor{s} = pcor_; errs{s} = errs_;
    clear Yp cors_ pcor_ errs_
end

%% visualize eval 
H = 2; W = length(Ypred); % [corr, err] * sys 
fig1 = figure('Units','normalized', 'Position',[.1,.1,.8,.8]);
for s = 1:W
    cors_ = cors{s}; errs_ = errs{s}; pcor_ = pcor{s};

    % plot correlation (top)
    subplot(H,W, s)
    t = (1:height(cors_)) - 1; % hzn (samples)
    t = seconds(t/fs); 
    plot(t, cors_(:,1:(end-1)))
    hold on; grid on; 
    plot(t, cors_(:,end), 'k', 'LineWidth',2) % aggregate 
    plot(t, mean(cors_(:,1:(end-1)),2), '--k', 'LineWidth',2) % avg
    title(sysName{s}); xlabel('Prediction Horizon');
    ylabel('Pearsons \rho')
    
    % plot error (bottom)
    subplot(H,W, W+s)
    t = (1:height(errs_)) - 1; % hzn (samples)
    t = seconds(t/fs); 
    plot(t, 100*errs_(:,1:(end-1)))
    hold on; grid on; 
    plot(t, 100*errs_(:,end), 'k', 'LineWidth',2) % aggregate 
    plot(t, 100*mean(errs_(:,1:(end-1)),2), '--k', 'LineWidth',2) % avg
    title(sysName{s}); xlabel('Prediction Horizon');
    ylabel('% RMSE')

    clear t cors_ errs_ pcor_
end

%% save 
svname = inputdlg('Save systems as:', 'File Save Name', 1, ...
    {[fn,'_simeval']});
if ~isempty(svname)
    svname = svname{1};
    save(fullfile(fp,[svname,'.mat']), 'sys', 'sysName', ...
        'dataTrain', 'dataTest', 'fn', ...'Ypred', ...
        'cors', 'pcor', 'errs')
    saveas(fig1, fullfile(fp,svname),'fig'); 
    saveas(fig1, fullfile(fp,svname),'png');
end