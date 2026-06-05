%% load data 
load('/Users/torenarginteanu/Desktop/Data_PD/PD26N003/Neuro Omega/SavedTable1375HzRT.mat')
Tbl = sortrows(Tbl, 'Time');
t = seconds(Tbl.Time);
x = Tbl.CLFP_AP_T___Central;
Fs = 1375; % Hz
load("/Users/torenarginteanu/Desktop/Data_PD/PD26N003/Neuro Omega/SPK_RT_SelectedTimes.mat")
load('/Users/torenarginteanu/Desktop/Data_PD/PD26N003/Neuro Omega/spikesort_227669d9c8fe69f115266d7020ddc4244c5f3779(2).mat')
spkTbl = depth_p0496_1;
tSpk = spkIdx/fs + seconds(spkTbl.Time(1));
ku = unique(kidx); tSpkK = cell(size(ku));
for ki = 1:length(ku)
    tSpkK{ki} = spkIdx(kidx == ku(ki)) / fs + seconds(spkTbl.Time(1));
end
tsel = (t >= seconds(spkTbl.Time(1))) & (t <= seconds(spkTbl.Time(end)));
x = x(tsel); t = t(tsel);

%% filter LFP 
hpf = fir1(1024,0.25/(Fs/2),"high");
xh = filtfilt(hpf,1,x); xl = x-xh;

%% build spike train 
z = zeros(size(t)); 
for zi = tSpk'
    [~,ti] = min(abs(t-zi));
    z(ti) = 1;
end
zk = cell(size(ku));
for ki = 1:length(ku)
    zk{ki} = zeros(size(t)); % Initialize spike train for each unique index
    for zi = tSpkK{ki}'
        [~, ti] = min(abs(t - zi));
        zk{ki}(ti) = 1; % Mark spikes in the corresponding train
    end
end

%% gaussian smooth pulse train 
%{
wvals = 1:ceil(100*Fs);
R = nan(size(wvals));
for wi = 1:length(wvals)
    w = wvals(wi);
    zw = smoothdata(z,1,'gaussian',w);
    R(wi) = corr(x, zw);
end
[~,wi] = max(R);
w = wvals(wi);
%}
w = 1*Fs;
zw = smoothdata(z,1,'gaussian',w);
zkw = cellfun(@(zi) smoothdata(zi,1,'gaussian',w), zk, 'UniformOutput',false);

%% compare spike and LFP signals 

% normalize 
xn = (x-mean(x))/std(x);
zkwn = cellfun(@(zi) (zi-mean(zi))/std(zi), zkw, 'UniformOutput',false);
zwn = (zw-mean(zw))/std(zw);

% evaluate 
r = corr(xn, zwn); 
rr = cellfun(@(zi) corr(xn,zi), zkwn);
lgd = string(ku)+": \rho="+string(rr);
lgd = ["raw LFP"; "all: \rho="+string(r); lgd];

% time domain 
figure; 
plot(t, xn); hold on; grid on; plot(t, zwn); 
for ki = 1:length(zkwn)
    plot(t, zkwn{ki}); % Plot each spike train for comparison
end
title('Spike-LFP time domain comparison');
xlabel('time (s)'); ylabel('normalized');
legend(lgd);

% freq domain 
figure; 
pwelch(xn,[],[],[],Fs,'power'); hold on;
pwelch(zwn,[],[],[],Fs,'power');
for ki = 1:length(zkwn)
    pwelch(zkwn{ki},[],[],[],Fs,'power');
end

%% modulated pulse train analysis 
fc = length(spkIdx)/(t(end)-t(1)); % initial est carrier freq = avg spk rate 
zw = movsum(z,Fs); % est inst freq

% compute Fourier series approximation of zw as sum of sine waves (with phase shifts)
% treat zw as real-valued signal sampled at Fs over times t
y = zw(:); Nt = numel(y);
Tdur = (t(end)-t(1)); % duration in seconds (duration object -> seconds)
dt = 1/Fs;
% prepare time vector in seconds
tv = (0:Nt-1)'*dt;

% compute FFT and corresponding freqs
Y = fft(y);
freqs = (0:Nt-1)'*(Fs/Nt);

% keep positive frequencies up to Nyquist
nPos = floor(Nt/2);
posIdx = 2:nPos; % exclude DC for sin representation (DC handled separately)
% amplitudes and phases for positive freqs
Apos = 2*abs(Y(posIdx))/Nt;        % amplitude of sine+cos pair
phpos = angle(Y(posIdx));         % phase of complex exponential

% include DC term separately
A0 = real(Y(1))/Nt; % DC (constant) term

% build sine-wave-only representation: A_k * sin(2*pi*f_k*t + phi_k)
fk = freqs(posIdx);
Ak = Apos;
phk = phpos - pi/2; % convert complex exponential phase to sine-phase
% (since e^{i(ωt+ph)} = cos(ωt+ph) + i sin(...); to express as sin(ωt+phi_s):
% sin(ωt+phi_s) = cos(ωt+phi_s - pi/2), adjust accordingly. Using phase shift:
% amplitude already accounts for factor 2 above.)

% reconstruct using a limited number of harmonics for stability (optional)
% choose K components with largest amplitudes
K = min(50, numel(Ak)); % limit to 50 components or less
[~,ord] = sort(Ak,'descend');
sel = sort(ord(1:K));
Ak = -Ak;

% build approximation
y_approx = A0 + zeros(size(tv));
for ii = 1:numel(sel)
    k = sel(ii);
    y_approx = y_approx + Ak(k)*sin(2*pi*fk(k)*tv + phk(k));
end

% provide outputs in workspace: fk(sel), Ak(sel), phk(sel), y_approx
fsine = fk(sel);
ampsine = Ak(sel);
phsine = phk(sel);

% optional plot to compare original and approximation
figure;
subplot(3,1,1);
plot(t, y); hold on; plot(t, y_approx); grid on;
xlabel('time (s)'); ylabel('inst freq (Hz)'); 
legend('original zw','sine-series approx');
title('Sine-wave Fourier series approximation (phase-shifted sines)');
subplot(3,1,2);
plot(freqs, 2*abs(Y)/Nt); hold on; stem(fsine, -ampsine);
ylabel('Amplitude (unitless)');
grid on; xlabel('frequency (Hz)'); 
xlim([min(fsine), max(fsine)])
subplot(3,1,3); 
plot(freqs, angle(Y)); hold on; stem(fsine, phsine + pi/2);
ylabel('Phase (rad)');
grid on; xlabel('frequency (Hz)'); 
xlim([min(fsine), max(fsine)])