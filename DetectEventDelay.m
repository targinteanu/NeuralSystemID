function [estDelay, EPs, fig] = DetectEventDelay(T, eventwindow, Tbase)
% Determine whether the timing events in table <T> are delayed by examining
% the evoked potentials within ± <eventwindow> seconds surrounding reported
% event times. If baseline values <Tbase> are provided, evoked potentials
% will be adjusted to reflect change from baseline.

% handle input args 
if nargin < 3
    Tbase = []; % bypass baseline adjustment
    if nargin < 2
        eventwindow = 0.5; % Default event window if not provided
    end
end
% output default args if unable to compute 
estDelay = nan; 
EPs = [];
fig = [];

% extract T properties 
if ~isempty(T)
fs = T.Properties.SampleRate;
if isnan(fs)
    fs = 1/median(seconds(diff(T.Time)));
end
L = floor(eventwindow*fs); % samples before and after 
Ev = T.Properties.Events;
if ~isempty(Ev)

% select events to be analyzed 
EvTypeSel = ...
    contains(lower(Ev.EventLabels), 'stim') | ...
    contains(lower(Ev.EventLabels), 'trig');
Ev = Ev(EvTypeSel, :);
Ev = Ev(Ev.Time > Ev.Time(1)+1.15*seconds(eventwindow), :);
Ev = Ev(Ev.Time < Ev.Time(end)-1.15*seconds(eventwindow), :);
if ~isempty(Ev)

% gather all evoked potentials 
EPs = nan(2*L+1, width(T), height(Ev));
for Ei = 1:height(Ev)
    [~,ti] = min(abs(T.Time - Ev.Time(Ei)));
    EPs(:,:,Ei) = T{ti+(-L:L),:};
end

% compute stats 
EPt = (-L:L)/fs; EPt = EPt'; % seconds 
EPavg = mean(EPs,3, 'omitnan');
EPstd = std(EPs,[],3, 'omitnan');
myttl = 'Signal ± 1SD';
if ~isempty(Tbase)
    baseAvg = mean(Tbase.Variables, 'omitnan');
    baseStd = std(Tbase.Variables, 'omitnan');
    EPavg = EPavg - baseAvg;
    EPstd = EPstd./baseStd;
    myttl = [myttl,', baseline-adjusted'];
end

% prepare figure 
EPpatch = [EPavg+EPstd; flipud(EPavg)-flipud(EPstd)];
fig = figure; 
for ch = 1:width(T)
    patch([EPt; flipud(EPt)], ...
        EPpatch(:,ch), ...
        colorwheel(ch/width(T)), 'FaceAlpha',.5, 'EdgeColor','none');
    hold on;
end
grid on;
xlabel('time from stim (s)'); 
ylabel(myttl)
legend(T.Properties.VariableNames, 'Location','westoutside');
title(T.Properties.Description);

end
end
end

end