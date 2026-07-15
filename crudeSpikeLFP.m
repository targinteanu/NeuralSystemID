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

[px,frange] = pwelch(xn,[],[],[],Fs,'power');

fc = length(spkIdx)/(t(end)-t(1)); % initial est carrier freq = avg spk rate 
zw = movsum(z,Fs); % est inst freq
zw0 = mean(zw); zw = zw-zw0; % different from carrier freq?
r = zw0/fc; % "repetition rate"

% inst freq fourier sine coeffs 
[fsine, ampsine, phsine, amp0] = FourierSine(zw, Fs, t);

% pulse train fourier coeffs 
wf = mean(WF); %wf = wf-mean(wf); % 0 offset
Tlen = ceil(fs/fc); % samples 
wfT = zeros(1,Tlen); wfT(1:length(wf)) = wf; % single period 
wft = 1:length(wfT); wft = (wft-1)/fs;
C0 = fc*sum(wfT)*(1/fs);
Ck = @(k) fc*sum(wfT.*exp(-1i*k*2*pi*fc*wft))*(1/fs);

% k=0 component 
%{
f0 = 2*pi*fsine'; f0 = [-fliplr(f0),0,f0];
q0 = [flipud(ampsine).*exp(-1i*flipud(phsine)); fc; ampsine.*exp(1i*phsine)]';
%}
f0 = [0;  fsine]';
q0 = [fc; ampsine.*exp(1i*phsine)/(2*r)]';

figure; 
plot(frange, 20*log10(px)); hold on; grid on; 
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
plot([f0,f1,f1b],20*log10(abs([q0,q1,q1b])),'.');

%% helpers

function [fk,qk] = getFQ(k, r, fc, fsine, ampsine, phsine)
% for k > 0 only! 

% bessel 
zrange = k*ampsine./(r*fsine); N = ceil(max(abs(zrange)));
nrange = -N:N;
J = zeros(length(nrange),length(zrange));
for ni = 1:length(nrange)
    J(ni,:) = besselj(nrange(ni),zrange); % for k=1
end
%{
% plot heatmap of Bessel matrix J vs nrange
figure;
imagesc(ampsine, nrange, J); % x: sorted amplitudes, y: nrange
axis xy;
colormap(parula);
colorbar;
xlabel('sine amplitudes');
ylabel('n (Bessel order)');
title('Bessel J_n(ampsine) vs n and amplitude');
%}
J(J.^2 < sum(J(:).^2)/numel(J)) = nan;

% first term (freq unshifted)
phshift = exp(1i*phsine*nrange)'; %phshiftAnti = exp(1i*(phsine+pi)*nrange)';
phsh2 = k*ampsine.*sin(phsine)./(r*fsine); phsh2 = sum(phsh2);
%phsh2Anti = exp(-1i*phsh2); phsh2 = exp(1i*phsh2);
Q1 = J.*phshift*phsh2;
numpk = prod(sum(~isnan(J)));
q1 = single(Q1(:,1)); f1 = single(fsine(1)*nrange');
f1 = f1(~isnan(q1))'; q1 = q1(~isnan(q1))';
if isempty(q1)
    f1 = 0; q1 = 1;
end
for l = 2:length(ampsine)
    ql = single(Q1(:,l)); fl = single(fsine(l)*nrange');
    fl = fl(~isnan(ql))'; ql = ql(~isnan(ql))';
    ff1 = (fl'+f1); f1 = single(ff1(:))';
    qq1 = ql'*q1; q1 = single(qq1(:))';
end
f1=f1+k*fc; q1 = q1*fc;

% second term (freq shifted)
fshift = [-flipud(fsine);fsine]';
Q1b = [flipud(ampsine).*exp(-1i*flipud(phsine)); ampsine.*exp(1i*phsine)]';
Q1b = reshape(Q1b,[1,1,length(Q1b)]); fshift = reshape(fshift,[1,1,length(fshift)]);
Q1b = Q1b.*Q1; F1b = (fsine*nrange)'.*fshift;
q1b = single(Q1b(:,1,:)); q1b = q1b(:);
f1b = single(F1b(:,1,:)); f1b = f1b(:);
f1b = f1b(~isnan(q1b))'; q1b = q1b(~isnan(q1b))';
if isempty(q1b)
    f1b = 0; q1b = 1;
end
for l = 2:length(ampsine)
    ql = single(Q1b(:,l,:)); fl = single(F1b(:,1,:)); ql = ql(:); fl = fl(:);
    fl = fl(~isnan(ql)); ql = ql(~isnan(ql))';
    if ~isempty(ql)
        ff1 = (fl+f1b); f1b = single(ff1(:))';
        qq1 = ql'*q1b; q1b = single(qq1(:))';
    end
end
f1b=f1b+k*fc; q1b = q1b/(2*r);

% output 
fk = [f1, f1b]; qk = [q1, q1b];

end


function [fsine, ampsine, phsine, A0] = FourierSine(zw, Fs, t)

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
K = min(8, numel(Ak)); % limit to 50 components or less
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

end