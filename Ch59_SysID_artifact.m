%% load, define 

load("Ch59.mat"); 
ind_bl_str = 4.5e5; ind_bl_end = 9.5e5; % no-stim baseline 
ind_rec_end = 2e6; % exclude dbs stim that is not 

Tr=NS2.Data(64,:);
Tr = Tr > 1e4;
Tr = Tr(1:ind_rec_end);
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

%% sys ID 

tic
[predBL, predAll, trnEval, tstEval, A] = fitLTIauton(dtaBL,dta);
toc
tstEval

KA = .0000001*eye(width(dta));
tic
[adaptBL,adaptAll,adaptTrnEval,adaptTstEval] = ...
    AID_LTI_auton(dtaBL,dta,[],KA,Tr);
toc
%{
tic
[adaptBL,adaptAll,adaptTrnEval,adaptTstEval] = ...
    AID_LTI_auton([],dtaBL,[],KA);
toc
%}
adaptTstEval

%% plot chan 59
chtoplot = dta.Properties.VariableNames{59};
figure; plot(dta, chtoplot, 'LineWidth',1); 
hold on; grid on; 
plot(predAll, chtoplot); plot(adaptAll, chtoplot);
legend('orig', 'fit', 'adapt')