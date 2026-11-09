function [rhos,Inds] = modCFCwindowed(x, Fs, bndMod, bndCar, winlen, winstride)

if nargin < 6
    winstride = [];
end
if isempty(winstride)
    winstride = 1;
end
L = length(x);
if winstride < 1
    error('Window stride must be positive.')
end
if winlen < 1
    error('window length must be positive.')
end

FMod = buildFIRBPF(Fs,bndMod(1),bndMod(2),2,1023);
FCar = buildFIRBPF(Fs,bndCar(1),bndCar(2),2,1023);
xMod = filtfilt(FMod,1,x);
xCar = filtfilt(FCar,1,x);

amp = envelope(xCar);
ampsmooth = filtfilt(FMod,1,amp);
%ampsmooth = movmean(amp, Fs/bndMod(2));

t0 = 1; tf = L;
% determine window start/ends
t1 = t0 + ceil(winlen/2); % define first center
Inds = t1:winstride:tf; % centers 
t2s = Inds + ceil(winlen/2);
Inds = Inds(t2s <= tf); t2s = t2s(t2s <= tf);
t1s = t2s - floor(winlen);
if isempty(Inds)
    error('signal is not long enough for window size.')
end

% perform PAC on windowed data 
rhos = nan(size(Inds));
for i = 1:length(Inds)
    wini = t1s(i):t2s(i);
    %rhos(i) = corr(xMod(wini), ampsmooth(wini));
    rr = xcorr(xMod(wini), ampsmooth(wini));
    rhos(i) = max(rr);
end

end