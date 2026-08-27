%% load raw data 
load('/Users/torenarginteanu/Desktop/Data_PD/PD26N002/Neuro Omega/SavedTable1375HzRT.mat')
Tbl = sortrows(Tbl, 'Time');
%x = Tbl.CLFP_AP_T___Central;
Fs = 1375; % Hz
load("/Users/torenarginteanu/Desktop/Data_PD/PD26N002/Neuro Omega/SPK_RT_SelectedTimes.mat")
spkTbl = TblD;
for xch = 1:width(Tbl)
    try
t = seconds(Tbl.Time);
tsel = (t >= seconds(spkTbl.Time(1))) & (t <= seconds(spkTbl.Time(end)));
x = Tbl{:,xch}; xname = Tbl.Properties.VariableNames{xch};
x = x(tsel); t = t(tsel);

% spectrogram; look for beta bursts 
figure; 
ax(1) = subplot(211);
spectrogram(x,5*Fs,[],[],Fs,"yaxis","power"); ylim([0 200]);
title(xname);
ax(2) = subplot(212);
plot(minutes(t), x); grid on; 
linkaxes(ax, 'x');
    catch ME
        warning(ME.message)
    end
end