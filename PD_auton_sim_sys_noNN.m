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
fsNew = dataTrain.Properties.SampleRate;

%% simulation 
N = width(sysAR.A{1,1}) - 1; % order of AR model
Cinv = pinv(bgLTI.C); % quick convert outputs to estimate of hidden states 
simSkipNum = 100; % new sim starts at every _th sample (for time/memory)
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

            % LTI simulations 
            sN = ceil(height(dataTest_)/simSkipNum);
            YsLTI = nan(height(dataTest_), width(dataTest_), sN );
            si = 1;
            for ti = 1:simSkipNum:height(dataTest_)
                disp(['LTI: simulation ',num2str(si),' of ',num2str(sN)])
                x0 = Cinv * dataTest_{ti,:}';
                ysLTI = sim(bgLTI, t, simOptions('InitialCondition',x0));
                YsLTI(:,:,si) = ysLTI.Variables;
                %t.Time = t.Time + t.Properties.TimeStep; % shift to next timestep
                si = si+1;
                clear ysLTI
            end

            % AR simulations 
            sN = ceil( (height(dataTest_)-N)/simSkipNum );
            YsAR = nan(height(dataTest_)-N+1, width(dataTest_), sN );
            si = 1;
            for ti = (N+1):simSkipNum:height(dataTest_)
                disp(['AR: simulation ',num2str(si),' of ',num2str(sN)])
                ysAR = myFastForecastAR(sysAR, dataTest_{(ti-N):(ti-1),:}, ...
                    height(dataTest_)-N);
                ysAR = [dataTest_{ti-1,:}; ysAR]; % ??
                YsAR(:,:,si) = ysAR;
                si = si+1;
                clear ysAR
            end

            Ysim{d,1} = YsLTI; Ysim{d,2} = YsAR;
            %clear YsLTI YsAR
        end
    end
    clear dataTest_
end

%% sim -> pred 
Ypred = cell(height(dataTest), width(sys)); % test set * sys * hzn

for d = 1:height(dataTest)
    dataTest_ = dataTest{d};
    if ~isempty(dataTest_)
        YsLTI = Ysim{d,1}; YsAR = Ysim{d,2};

        YpLTI = cell(height(YsLTI),1);
        for h = 1:height(YsLTI)
            y = dataTest_{h:simSkipNum:end, :}; 
            y = y(1:ceil(height(YsLTI)/(simSkipNum)), :); 
            yp = squeeze(YsLTI(h,:,:))';
            YpLTI{h} = cat(3,y,yp);
            clear y yp
        end

        YpAR = cell(height(YsAR),1);
        for h = 1:height(YsAR)
            y = dataTest_{h:simSkipNum:end, :}; 
            y = y(1:ceil(height(YsAR)/(simSkipNum)), :); 
            %y = y(N:end,:);
            yp = squeeze(YsAR(h,:,:))';
            YpAR{h} = cat(3,y,yp);
            clear y yp
        end

        Ypred{d,1} = YpLTI; Ypred{d,2} = YpAR;
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

%% evaluate pred 
Yp = Ypred{1};
cors = nan(height(Yp), width(dataTrain)+1); % hzn * chan (incl overall)
pcor = cors;
for h = 1:height(cors)
    for c = 1:width(dataTrain)
        y = Yp{h}(:,c,1); yp = Yp{h}(:,c,2);
        [cors(h,c), pcor(h,c)] = corr(y,yp,"tail","right");
    end
    y = Yp{h}(:,:,1); yp = Yp{h}(:,:,2);
    y = y(:); yp = yp(:);
    [cors(h,c+1), pcor(h,c+1)] = corr(y,yp,"tail","right");
    clear y yp
end