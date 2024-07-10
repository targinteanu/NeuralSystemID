%% neural sys ID - train
%{
dta_nssTrain = dtaBL;
dta_nssTrian = resample(dta_nssTrain, 1, 4); % downsample
dta_nssTrain = dta_nssTrian(50000:100000,:);
nss_str = idNeuralStateSpace(width(dta));
nss_str.OutputName = dtaBL.Properties.VariableNames;
nss_trn = nlssest(dta_nssTrain,nss_str);
%}
load('NeuralSysID.mat');
%% neural sys ID - evaluate
chtoplot = dta.Properties.VariableNames{59};

[x,y,xLin,yLin] = evaluateTimeTableAuton(nss_trn, dtaBL);

% training prediction 
opt = simOptions('OutputTimes',seconds(dta_nssTrain.Time));
dtaBLest = sim(nss_trn,dta_nssTrain(1,:),opt);
figure; plot(dta_nssTrain, chtoplot); grid on; hold on; 
plot(dtaBLest, chtoplot);

% training prediction 
dtaBLest = evaluate(nss_trn,table2array(dta_nssTrain));
figure; plot(dta_nssTrain, chtoplot); grid on; hold on; 
plot(dtaBLest, chtoplot);

% baseline prediction 
opt = simOptions('OutputTimes',seconds(dtaBL.Time));
dtaBLest = sim(nss_trn,dtaBL(1,:),opt);
figure; plot(dtaBL, chtoplot); grid on; hold on; 
plot(dtaBLest, chtoplot);

% complete prediction 
opt = simOptions('OutputTimes',seconds(dta.Time));
dta_est = sim(nss_trn,dta(1,:),opt);
figure; plot(dta, chtoplot); grid on; hold on; 
plot(dta_est, chtoplot);