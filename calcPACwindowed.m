function [MIs,Inds] = calcPACwindowed(xPhi, xAmp, winlen, winstride, nbin)
% 
% Get phase-amplitude-coupling (PAC) modulation index (MI) in rolling
% windows. 
% 
% Inputs: 
%   xPhi: phase signal in time domain, filtered 
%   xAmp: amp signal in time domain, filtered
%   winlen: length of window for PAC calculations (in samples) 
%   winstride: distance between windows, i.e. 1-overlap, default 1 (samples)
%   nbin: # of phase bins, default 18
% 
% Outputs: 
%   MIs: modulation index (MI) in windows 
%   Inds: sample indexes (of orig. signals) corresponding to MIs 
% 

if nargin < 5
    nbin = [];
    if nargin < 4
        winstride = [];
    end
end
if isempty(nbin)
    nbin = 18;
end
if isempty(winstride)
    winstride = 1;
end
L = length(xPhi);
if L ~= length(xAmp)
    error('Signals need to be same length and aligned in time.')
end
if nbin < 1
    error('Wrong number of bins requested.')
end
if winstride < 1
    error('Window stride must be positive.')
end
if winlen < 1
    error('window length must be positive.')
end

% hilbert 
phi = angle(hilbert(xPhi));
amp = abs(hilbert(xAmp));

% focus on the central 80% of data to avoid edge effects
t0 = floor(.1*L); tf = ceil(.9*L);

% determine window start/ends
t1 = t0 + ceil(winlen/2); % define first center
Inds = t1:winstride:tf; % centers 
t2s = Inds + ceil(winlen/2);
Inds = Inds(t2s <= tf); t2s = t2s(t2s <= tf);
t1s = Inds - floor(winlen);
if isempty(Inds)
    error('signal is not long enough for window size.')
end

% perform PAC on windowed data 
MIs = nan(size(Inds));
for i = 1:length(Inds)
    wini = t1s(i):t2s(i);
    MIs(i) = calcPAChelper(phi(wini), amp(wini), nbin, false);
    if isnan(MIs(i)) || isinf(MIs(i))
        %keyboard
    end
end

end