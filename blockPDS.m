function [t2phi, i2phi, phi_inst, f_inst, A_avg] = ...
    blockPDS(pastData, futureData, fs, phi, tmin, fmin, fmax, usefuture)
% Determine the time (s) and # samples to next desired phase phi from a
% block of data sampled at a constant rate fs (Hz). Also return the current
% inst. phase phi_inst (rad) and frequency f_inst (Hz).
% Block data should include some length of pastData and
% (forecasted/predicted) futureData to minimize edge effects at the present
% timepoint, which is the last element of pastData. Data should be input as
% columns. 
% phi [desired] is in radians, i.e. phi=0 for peak, phi=pi for trough
% frequency will be clipped within range [fmin, fmax] (Hz) 
% tmin is the minimum time ahead (s) to the desired phase; i.e. ignore
% targets closer to this that are within the expected system delay
% usefuture flag: if true, futureData will be will be hilbert transformed
% and taken as truth to identify target phase; otherwise, target will be
% estimated uding instantaneous phase and freq. and assuming constant freq.

if nargin < 8
    usefuture = true;
    if nargin < 7
        fmax = inf;
        if nargin < 6
            fmin = -inf;
            if nargin < 5
                tmin = 0;
            end
        end
    end
end

N = size(pastData,1); M = size(futureData,1);
blockData = [pastData; futureData];

[phi_block, f_block, A_block] = instPhaseFreq(blockData, fs);
A_avg = mean(A_block);
phi_inst = phi_block(N,:);
f_block = max(f_block, fmin); 
f_block = min(f_block, fmax);
fwinlen = floor(.1*N); 
%fwinlen = min(fwinlen, M);
%fwin = N + ((-fwinlen):fwinlen);
%f_inst = mean(f_block(fwin,:));
f_block = smoothdata(f_block, "gaussian", fwinlen);
f_inst = f_block(N);
%[~,f_inst] = zerocrossrate(blockData); f_inst = f_inst*fs/(2*length(blockData));
T=1/f_inst;

% examine the next full phase cycle
imin = ceil(tmin*fs); 
if usefuture
phi_future = phi_block((N+imin):end,:);
phi_fut_1 = phi_future(1); phi_fut_ = unwrap(phi_future);
[~,phi_fut_end] = min(abs(phi_fut_-(2*pi + phi_fut_1)));
phi_fut_end = min(length(phi_future), phi_fut_end+1);
phi_future = phi_future(1:phi_fut_end);
end

% time to next [desired] phi 
t2phi = zeros(size(phi)); i2phi = t2phi;
for p = 1:length(phi)
    phi_ = phi(p);
    if isnan(phi_)
        t2phi(p) = inf;
        i2phi(p) = inf;
    else

    if usefuture
    % accurate & full-cycle AR model assumption 
    [~,i] = min(abs(radfix(phi_future-phi_)));
    i = i + imin;
    i2phi(p) = i; 
    t2phi(p) = i/fs;

    else
    % const freq assumption
    t = (mod(phi_+2*pi-phi_inst,2*pi)./f_inst)/(2*pi); 
    % account for minimum delay time tmin 
    if t < tmin
        nT = (tmin-t)/T; % how many periods needed to add 
        t = t + ceil(nT)*T; 
    end
    t2phi(p) = t;
    i2phi(p) = floor(fs*t2phi(p));
    
    end
    end
end