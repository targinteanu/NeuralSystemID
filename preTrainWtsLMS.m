function [w, e, op, fig] = preTrainWtsLMS(g, d, N, nUpdates, dLMS)
% Use a subset of the signal to train optimal LMS weights. 
%
% Inputs: 
%   g: training noise reference signal as single column corresponding to d 
%   d: training unfiltered signal in timetable 
%   N: number of filter taps 
%   nUpdates: how many times to display progress and whether or not to plot weights. 
%             0 = no output and no plots 
%             Default = 10
%   dLMS: if false, does LMS; if true, does dLMS. Default false
% 
% Outputs: 
%   w: trained weights, as columns, from oldest --> current time 
%   e: LMS error signal for training period 
%   op: LMS output signal for training period 
%   fig: matlab figure with the trained weights of each channel

if nargin < 5
    dLMS = false;
end

% TO DO: account for time discontinuity 

dlen = height(d); numchan = width(d);
chname = d.Properties.VariableNames;
t = d.Time;
e = nan(dlen-N+1, numchan);
t_e = t(N:end);
d = table2array(d);

if dlen ~= length(g)
    error('d and g must be the same length and duration.')
end

%% organize training epochs 
G = zeros(dlen-N+1, N, numchan); 
D = zeros(dlen-N+1, numchan);
for idx = 1:numchan
    ch = chname{idx};
    D(:,idx) = d(N:dlen, idx);
    for nf = 1:(dlen-N+1)
        if nUpdates
            if ~mod(nf, floor(dlen/(nUpdates)))
                disp(['Building Channel ',ch,' Training Matrix: ',num2str(100*nf/dlen),'%']);
            end
        end
        G(nf,:,idx) = g(nf:(nf+N-1));
    end
end

%% training  
if nUpdates
    fig = figure('Units','normalized', 'Position',[.1 .1 .4 .8]);
end
w = zeros(N, numchan);
for idx = 1:numchan
    ch = chname{idx};
    Gidx = G(:,:,idx); Didx = D(:,idx);
    dG = diff(Gidx); dD = diff(Didx);
    if ~dLMS
        w(:,idx) = (((Gidx'*Gidx)^-1)*Gidx')*Didx;
    else
        w(:,idx) = (((dG'*dG)^-1)*dG')*dD;
    end
    if nUpdates
        figure(fig); subplot(numchan, 1, idx); stem(w(:,idx)); grid on;
        title(['Channel ',ch,' training']);
        xlabel('tap'); ylabel('weight');
        pause(eps);
    end
end
if nUpdates
    pause(.5);
end

wnan = isnan(w);
if sum(wnan(:))
    w(wnan) = 0;
    wnan = sum(wnan(:))/numel(wnan);
    warning(['Weights are ',num2str(wnan*100),'% undefined; reverting to zeros.'])
end

%% post-processing  
op = zeros([dlen-N+1,numchan]); 
for idx = 1:numchan
    op(:,idx) = G(:,:,idx)     *w(:,idx);
end

e = d; e(N:end,:) = e(N:end,:) - op;

e  = array2timetable(e,  "RowTimes",t,   "VariableNames",chname);
op = array2timetable(op, "RowTimes",t_e, "VariableNames",chname);

end