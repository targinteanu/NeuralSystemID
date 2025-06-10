%% load, define 

load("Ch59.mat"); 
ind_bl_str = 4.5e5; ind_bl_end = 9.5e5; % no-stim baseline 
ind_rec_end = 2e6; % exclude dbs stim that is not annotated

Tr=NS2.Data(64,:);
Tr = Tr(1:ind_rec_end);
Tr_thr = Tr > 1e4;
nDelay = 3;
Tr_thr = [Tr_thr((nDelay+1):end), false(1,nDelay)]; % start nDelay samples earlier
Fs = NS2.MetaTags.SamplingFreq;

dta = ns2timetable(NS2); 
dta = dta(:,1:63); % exclude stim (analog in)
dta = dta(:,[1:2,4:53,55:56,58:60,62:end]); % exclude stim/noisy/ref channels from analysis
dta = dta(1:ind_rec_end,:); 

chnum_to_plot = 56; % should be LS59
chtoplot = dta.Properties.VariableNames{chnum_to_plot}

%% preprocess 

%Fs = dta.Properties.SampleRate;
hpf = designfilt('highpassiir', 'SampleRate',Fs, 'DesignMethod','butter', ...
    'StopbandFrequency',.5, 'PassbandFrequency',1.5, ...
    'StopbandAttenuation',60, 'PassbandRipple',1);
dta = FilterTimetable(@(d,x) filtfilt(d,x), hpf, dta);
dtaBL = dta(ind_bl_str:ind_bl_end,:);

%% AR model 

dta1ch = dta(:,chnum_to_plot);
dtaBL1ch = dtaBL(:,chnum_to_plot);

% best fit and dynamic AR
N = 10; % order 
tic
[aarBL,aarAll,aarTrnEval,aarTstEval,ARmdl] = ...
    AdaptAR(dtaBL1ch,dta1ch,N,1e-12,Tr_thr,true);
[arSO,arEval] = projAR(ARmdl,dta1ch,Tr_thr);
toc
aarTstEval
arEval

%% LTI fit sys ID 

tic
[predBL, predAll, trnEval, tstEval, A] = fitLTIauton(dtaBL,dta);
[predSO,evalSO] = projLTIauton(A,dta,Tr_thr);
toc
tstEval
evalSO

% use training residual to estimate process noise Q
res = predBL - dtaBL; % residual 
Qcov = cov(res.Variables);

%% adaptive Sys ID

KA = (1e-4)*eye(width(dta));
Am = (-1e3)*eye(width(dta));
tic
[adaptBL,adaptAll,adaptTrnEval,adaptTstEval] = ...
    AID_LTI_auton(dtaBL,dta,Am,KA,Tr_thr,true);
toc
adaptTstEval

%% LMS filter artifact removal 
N = 20; % filter # taps

dta59 = dta(:,chnum_to_plot);
g = diff(double(Tr))'; g = g.*(g > 1e4); 
g = [g((nDelay+1):end); false(nDelay,1)]; % start nDelay samples earlier

% find the optimal weights 
ind_stim_str = find(Tr_thr); 
ind_stim_end = ind_stim_str(end); ind_stim_str = ind_stim_str(1);
ind_stim_str = max(1, ind_stim_str-N);
ind_stim_end = min(ind_rec_end, ind_stim_end+N);
g_noise_train = g(ind_stim_str:ind_stim_end);
dta_noise_train = dta59(ind_stim_str:ind_stim_end,:);
w0 = preTrainWtsLMS(g_noise_train,dta_noise_train,N,2,false);

dta59 = dta59(1000000:end,:);
g = g((1000000-1):end);
dtaLMS  = filterLMS(g,dta59,1,N,w0,100,false,true);
%dtaLMSd = filterLMS(double(Tr(1000000:end))',dta59,.01,32,[],100,true, true);

%% Kalman filter 

% characterize process noise waveform 
noise_ind = find(Tr_thr); 
noise_ind_new = diff(noise_ind) > 2;
noise_ind_new = find(noise_ind_new) + 1;
noise_inds = {};
for ind = 2:length(noise_ind_new)
    ind_curr = noise_ind_new(ind-1) : (noise_ind_new(ind)-1);
    noise_inds = [noise_inds, noise_ind(ind_curr)];
end
clear noise_ind ind ind_curr noise_ind_new
noise_tbls = cellfun(@(inds) dta(inds,:), noise_inds, 'UniformOutput',false);
figure; hold on; grid on;
for tbl = noise_tbls
    tbl = tbl{:};
    tbl.Time = tbl.Time - tbl.Time(1);
    plot(tbl,chtoplot, 'Color','r');
end
title('artifacts')

% evaluate measurement noise (co)variance 
N = 19; % samples 
Rt = tbl.Time(1:N);
RR = nan(width(dta),length(noise_tbls),N);
for m = 1:length(noise_tbls)
    tbl = noise_tbls{m};
    for n = 1:N
        RR(:,m,n) = tbl{n,:};
    end
end
Rstd = std(RR,[],2); Rstd = squeeze(Rstd);
Ravg = mean(RR,2); Ravg = squeeze(Ravg);
ravg = Ravg(chnum_to_plot,:);
plot(Rt,ravg, 'k', 'LineWidth',1);
plot(Rt,ravg+Rstd(chnum_to_plot,:), '--k', 'LineWidth',1);
plot(Rt,ravg-Rstd(chnum_to_plot,:), '--k', 'LineWidth',1);
Rcov = nan(width(dta),width(dta),N);
for n = 1:N
    Rcov(:,:,n) = cov(RR(:,:,n)');
end
Rvar = nan(size(Rstd)); 
for n = 1:N
    Rvar(:,n) = diag(Rcov(:,:,n));
end
plot(Rt,ravg'+sqrt(squeeze(Rcov(chnum_to_plot,chnum_to_plot,:))), ...
    ':b', 'LineWidth',1);
rmse(Rvar(:), Rstd(:).^2)

% run filter on data 
g = [0, diff(double(Tr))]'; g = g.*(g > 1e4);
g = [g((nDelay+1):end); false(nDelay,1)]; % start nDelay samples earlier
dtaKal = AdaptKalmanAuton(g,dta,Am,KA,A,1,Qcov,Rcov,[],true);

%% plot chan 59
figure; plot(dta, chtoplot, 'LineWidth',1); ybnd = ylim;
hold on; grid on; 
%plot(predAll, chtoplot); 
plot(predSO,chtoplot); 
plot(adaptAll, chtoplot);
plot(dtaLMS, chtoplot); 
%plot(dtaLMSd, chtoplot);
plot(dtaKal, chtoplot);
legend('orig', 'LTI', 'AID', 'LMS', 'Kal'); title(chtoplot);
ylim(ybnd); 