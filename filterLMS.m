function [e, w, fig] = filterLMS(g, d, stepsize, N, w, nUpdates, dLMS, nLMS)
% Perform online LMS adaptive filtering. 
% 
% Inputs: 
%   g: noise reference signal as single column corresponding to d
%   d: unfiltered signal in timetable 
%   stepsize: LMS step size / learning rate
%   N: number of filter taps 
%   w: starting weights 
%      Default = zeros
%   nUpdates: how many times to display progress and whether or not to plot weights 
%             and smoothed error. 0 = no output and no plots 
%             Default = 100
%   dLMS: if false, does LMS; if true, does dLMS. Default false
%   nLMS: if true, does NLMS; if false, does not normalize. Default false
% 
% Outputs: 
%   e: LMS error signal 
%   w: final weights 
%   fig: matlab figure showing the resulting weights and error signal 
 
if nargin < 8
    nLMS = false;
    if nargin < 7
        dLMS = false;
        if nargin < 6
            nUpdates = 100;
            if nargin < 5
                w = [];
            end
        end
    end
end
if isempty(w)
    w = zeros(N,width(d));
end

if nUpdates
    fig = figure('Units','normalized', 'Position',[.1 .1 .8 .8]);
end

dlen = height(d); numchan = width(d);
chname = d.Properties.VariableNames;
t = d.Time;
e = nan(dlen-N+1, numchan);
t_e = t(N:end);
d = table2array(d);

if dlen ~= length(g)
    error('d and g must be the same length and duration.')
end

for idx = 1:numchan
    ch = chname{idx};
    % train w: iterate grad descent
    if nUpdates
        figure(fig);
        subplot(numchan,2,2*idx-1); wplot = stem(w(:,idx));    grid on;
        title(['Channel ',ch,' online']);
        xlabel('tap'); ylabel('weight');
        subplot(numchan,2,2*idx);   eplot = semilogy(e(:,idx)); grid on;
        title(['Channel ',ch,' online']);
        xlabel('timepoint'); 
        if ~dLMS
            ylabel('e^2');
        else
            ylabel('(de/dt)^2');
        end
        pause(.5);
    end
    Gprev = zeros(1,N); Eprev = 0;
    for ep = (N:dlen)-N+1
        Gidx = g((1:N)+ep-1)';
        E = d(ep+N-1,idx) - Gidx*w(:,idx);
        e(ep, idx) = E;
        if ~dLMS
            dw = E*Gidx'; % LMS
            if nLMS
                dw = dw./(Gidx*Gidx' + eps);
            end
        else
            dG = Gidx-Gprev;
            dw = (E-Eprev) * dG'; % dLMS
            if nLMS
                dw = dw./(dG*dG' + eps);
            end
        end
        w(:,idx) = w(:,idx) + stepsize*dw;
        Gprev = Gidx; Eprev = E;
        if nUpdates
            if ~mod(ep, floor(dlen/nUpdates))
                wplot.YData = w(:,idx); 
                edata = e(:,idx);
                if dLMS
                    edata = [0;diff(edata)];
                end
                eplot.YData = movmean(edata.^2, 5000);
                disp(['Online Channel ',ch,': ',num2str(100*ep/dlen),'%'])
                pause(eps);
            end
        end
    end
end

e = array2timetable(e, "RowTimes",t_e, "VariableNames",chname);

end