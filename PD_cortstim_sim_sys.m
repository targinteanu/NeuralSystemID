%% Parkinson's Disease (PD) Project - simulate systems with cortical input
% Simulate at every time point in order to figure out how many samples
% ahead the prediction is valid. 

%% load data file 
[fn,fp] = uigetfile('*cortstim*.mat');
load(fullfile(fp,fn), 'bgLTI_A', 'bgLTI_NA', 'bgZIRZSR', 'bgHWnn', ...
    'dataTrain', 'dataTest');
disp([fp,' --- ',fn]);
[~,fn] = fileparts(fn);
sysName = {'bgLTI_A', 'bgLTI_NA', 'bgZIRZSR', 'bgHWnn'};
sys = cellfun(@eval,sysName, 'UniformOutput',false);
for s = 1:length(sysName)
    sysName{s}(sysName{s} == '_') = ' ';
end
fs = dataTrain.Properties.SampleRate;

%% simulation 
% Ysim will have entries for each system (AR and LTI).
% Entries have timepoints as rows and channels as
% columns. Simulations starting at regular intervals are stacked in the
% third dimension. 

simSkipNum = 10; % new sim starts at every _th sample (for time/memory)
simMaxDur = 5000; % samples; 75 sec
Ysim = cell(1, width(sys)); % all sim results 

% prepare starting point data for each set of test data 
dataTest_ = dataTest;
if ~isempty(dataTest_)
    dataTest_ = dataTest_(1:floor(height(dataTest_)/2), :);
    if ~isempty(dataTest_)
        dataTest_ = retime(dataTest_, 'regular','nearest',...
            'TimeStep', seconds(bgLTI_NA.Ts) );
        dataTest_.Time = dataTest_.Time - dataTest_.Time(1);
        u = dataTest_(:,end); % timetable of input only
        u = u(1:min(simMaxDur, height(u)), :);
        dataTest_ = dataTest_(:,1:(end-1));

        for s = 1:length(sysName)
            S = sys{s};
            sN = ceil(height(dataTest_)/simSkipNum);
            Ys = nan(height(u), width(dataTest_), sN );
            si = 1;
            for ti = 1:simSkipNum:height(dataTest_)
                if ~mod(si, 10)
                    disp([sysName{s},': simulation ',num2str(si),' of ',num2str(sN)])
                end
                if s < 4
                    x0 = pinv(S.C) * dataTest_{ti,:}';
                else
                    x0 = 'z';
                end
                ys = sim(S, u, simOptions('InitialCondition',x0));
                Ys(:,:,si) = ys.Variables;
                si = si+1;
                clear ys
            end

            Ysim{s} = Ys; 
            clear Ys S
        end

    end
end
clear dataTest_

%% sim -> pred 
% Ypred entries will have rows corresponding to prediction horizons, and
% these entries should be similar to a downsampled output of MATLAB's
% predict function ("yp") at that horizon, stacked along the third
% dimension with the corresponding actual data ("y").

Ypred = cell(1, width(sys)); % test set * sys * hzn

dataTest_ = dataTest;
if ~isempty(dataTest_)
    dataTest_ = dataTest_(:,1:(end-1)); % remove input 
    for s = 1:length(sysName)
        Ys = Ysim{s}; 

        Yp = cell(height(Ys),1);
        for h = 1:height(Ys)
            y = dataTest_{h:simSkipNum:end, :}; 
            yp = squeeze(Ys(h,:,:))';
            y = y(1:height(yp), :);
            Yp{h} = single(cat(3,y,yp));
            clear y yp
        end

        Ypred{s} = Yp; 
        clear Yp Ys S
    end
end

%% evaluate pred correlation and error
% In each channel (and across all channels together), evaluate the
% correlation and the fraction RMS error between the predicted and actual
% data as a function of prediction horizon. 

cors = cell(size(Ypred)); pcor = cors; % correlation and p-value 
errs = cell(size(Ypred)); % pRMSE 

for s = 1:length(Ypred)
    disp(['Evaluating ',sysName{s}])
    Yp = Ypred{s};
    cors_ = nan(height(Yp), width(dataTrain)); % hzn * chan (incl overall)
    pcor_ = cors_; errs_ = cors_;
    for h = 1:height(cors_) % hzn
        if ~mod(h, 50)
            disp([' - timepoint ',num2str(h),' of ',num2str(height(cors_))])
        end

        % eval each channel
        for c = 1:(width(dataTrain)-1)
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