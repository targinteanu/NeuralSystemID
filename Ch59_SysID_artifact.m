%% load, define 

load("Ch59.mat"); 
ind_bl_str = 4.5e5; ind_bl_end = 9.5e5; % no-stim baseline 
ind_rec_end = 2e6; % exclude dbs stim that is not 

Tr=NS2.Data(64,:);
Tr = Tr(1:ind_rec_end);
Tr_thr = Tr > 1e4;
Fs = NS2.MetaTags.SamplingFreq;

dta = ns2timetable(NS2); 
dta = dta(:,1:63); % exclude stim (analog in)
dta = dta(1:ind_rec_end,:); 

%% preprocess 

%Fs = dta.Properties.SampleRate;
hpf = designfilt('highpassiir', 'SampleRate',Fs, 'DesignMethod','butter', ...
    'StopbandFrequency',.5, 'PassbandFrequency',1.5, ...
    'StopbandAttenuation',60, 'PassbandRipple',1);
dta = FilterTimetable(@(d,x) filtfilt(d,x), hpf, dta);
dtaBL = dta(ind_bl_str:ind_bl_end,:);

%% LTI fit sys ID 

tic
[predBL, predAll, trnEval, tstEval, A] = fitLTIauton(dtaBL,dta);
[predSO,evalSO] = projLTIauton(A,dta,Tr_thr);
toc
tstEval
evalSO

%% adaptive Sys ID

KA = (1e-4)*eye(width(dta));
Am = (-1e3)*eye(width(dta));
tic
[adaptBL,adaptAll,adaptTrnEval,adaptTstEval] = ...
    AID_LTI_auton(dtaBL,dta,Am,KA,Tr_thr,false);
toc
adaptTstEval

%% neural sys ID 
nss_str = idNeuralStateSpace(width(dta));
nss_str.OutputName = dtaBL.Properties.VariableNames;
nss_trn = nlssest(dtaBL,nss_str);

%% LMS filter artifact removal 

dta59 = dta(:,59);
g = diff(double(Tr))'; g = g.*(g > 1e4); 

% find the optimal weights 
ind_stim_str = find(Tr_thr); 
ind_stim_end = ind_stim_str(end); ind_stim_str = ind_stim_str(1);
ind_stim_str = max(1, ind_stim_str-16);
ind_stim_end = min(ind_rec_end, ind_stim_end+16);
g_noise_train = g(ind_stim_str:ind_stim_end);
dta_noise_train = dta59(ind_stim_str:ind_stim_end,:);
w0 = preTrainWtsLMS(g_noise_train,dta_noise_train,16,2,false);

dta59 = dta59(1000000:end,:);
g = g((1000000-1):end);
dtaLMS  = filterLMS(g,dta59,1,16,w0,100,false,true);
%dtaLMSd = filterLMS(double(Tr(1000000:end))',dta59,.01,32,[],100,true, true);

%% plot chan 59
chtoplot = dta.Properties.VariableNames{59};
figure; plot(dta, chtoplot, 'LineWidth',1); ybnd = ylim;
hold on; grid on; 
%plot(predAll, chtoplot); 
%plot(predSO,chtoplot); 
plot(adaptAll, chtoplot);
plot(dtaLMS, chtoplot); 
%plot(dtaLMSd, chtoplot);
legend('orig', 'AID', 'LMS'); title(chtoplot);
ylim(ybnd); 