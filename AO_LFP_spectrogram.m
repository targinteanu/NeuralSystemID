%% load raw data 

load('/Users/torenarginteanu/Desktop/Data_PD/PD26N003/Neuro Omega/SavedTable1375HzLT.mat')
Tbl = sortrows(Tbl, 'Time');
%x = Tbl.CLFP_AP_T___Central;
Fs = 1375; % Hz
load("/Users/torenarginteanu/Desktop/Data_PD/PD26N003/Neuro Omega/SPK_LT_SelectedTimes(2).mat")
spkTbl = depth_p1886;

E = Tbl.Properties.Events;
E = E(contains(E.EventLabels, 'Stim'), :);
D = AlphaOmegaTable2Depth(Tbl);

% truncate time to match spk
t = (Tbl.Time);
tsel = (t >= (spkTbl.Time(1))) & (t <= (spkTbl.Time(end)));
t = t(tsel);
Tbl = Tbl(tsel,:); 
Esel = (E.Time >= (spkTbl.Time(1))) & (E.Time <= (spkTbl.Time(end)));
E = E(Esel,:);
Dsel = (D.Time >= (spkTbl.Time(1))) & (D.Time <= (spkTbl.Time(end)));
D = D(Dsel,:);

%%
for xch = 1:width(Tbl)
    try
x = Tbl{:,xch}; xname = Tbl.Properties.VariableNames{xch};

%% spectrogram; look for beta bursts 
figure; 
ax(1) = subplot(311);
spectrogram(x,1*Fs,[],[],Fs,"yaxis","power"); ylim([0 200]);
title(xname);
linkfirstax = true;
if contains(ax(1).XLabel.String, '(s)')
    timefun = @seconds;
elseif contains(ax(1).XLabel.String, '(minutes)')
    timefun = @minutes;
elseif contains(ax(1).XLabel.String, '(hours)')
    timefun = @hours;
else
    timefun = @(a) a; linkfirstax = false;
end
ax(2) = subplot(312);
plot(timefun(t-t(1)), x); grid on; ylabel('LFP sig')
ax(3) = subplot(313);
if ~isempty(D)
plot(timefun(D.Time-D.Time(1)), D.DEPTH); ylabel('Depth (mm)');
end
hold on; grid on; 
if ~isempty(E)
stem(timefun(E.Time-E.Time(1)), -5*ones(size(E)), '.'); 
end
if ~isempty(D) && ~isempty(E)
    legend('Depth', 'Stim');
end
if linkfirstax
    linkaxes(ax, 'x');
else
    linkaxes(ax(2:end), 'x');
end
    catch ME
        warning(ME.message)
    end
end