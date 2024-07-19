%% start EEGlab 
    eeglabpath = '/Applications/MATLAB_R2021b.app/toolbox/eeglab2022.0';
    addpath(eeglabpath);
    eeglab

%% load EEG 
basedir = cd; 
cd('/Users/torenarginteanu/Desktop/Data_Chronic Pain'); 
[fn,fp] = uigetfile('.mat', 'Select Pre- or Post-processed EEG');
cd(basedir); 
load([fp,filesep,fn]);

%% epoch 
[EEG0, Epo0] = epochBaseline(EEG_table,'BaselineOpen','before experiment',1,1);
[EEG1, Epo1] = epochStim(EEG_table,'PinPrick',1);

%% plot fit
plotchan = 'CZ'; 
A0 = plotModelFit(EEG0, Epo0, @(tsTbl, trTbl) fitLTIauton(tsTbl, trTbl), plotchan);
A1 = plotModelFit(EEG1, Epo1, @(tsTbl, trTbl) fitLTIauton(tsTbl, trTbl), plotchan);

%% plot source-sink
chanlocs = EEG1(1).chanlocs;
[srcness0,snkness0,figs] = SourceSink(A0, chanlocs);
figure(figs(2)); 
sgtitle({fn,'Baseline'});
[srcness1,snkness1,figs] = SourceSink(A1, chanlocs);
figure(figs(2)); 
sgtitle({fn,'Stim'});

% order chans as columns 
srcness0 = squeeze(srcness0)'; srcness1 = squeeze(srcness1)';
snkness0 = squeeze(snkness0)'; snkness1 = squeeze(snkness1)';

%% compare source-sink stim to baseline 
[Hsnk,Psnk,Csnk,Ssnk] = ttest2(snkness1,snkness0,"Vartype","unequal");
[Hsrc,Psrc,Csrc,Ssrc] = ttest2(srcness1,srcness0,"Vartype","unequal");
Tsnk = Ssnk.tstat; Tsrc = Ssrc.tstat; 
Trng = 'absmax'; % change to be same for source and sink?
figure; sgtitle(fn); 
subplot(2,2,1); topoplot(Tsrc, chanlocs, 'maplimits', Trng); 
title('Source T'); colorbar; 
subplot(2,2,2); topoplot(Tsnk, chanlocs, 'maplimits', Trng);
title('Sink T'); colorbar; 
subplot(2,2,3); topoplot(Psrc, chanlocs, 'maplimits', [0,1]);
title('Source p'); colorbar; 
subplot(2,2,4); topoplot(Psnk, chanlocs, 'maplimits', [0,1]);
title('Sink p'); colorbar;