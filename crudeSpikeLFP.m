%% load raw data 
minFreq = 0.3; % Hz
load('/Users/torenarginteanu/Desktop/Data_PD/PD24N009/Neuro Omega/mat/LT/Saved To Table/Table Data 2026-09-04 16.02.54 d0a6723833db66/SavedTable1375HzLT.mat')
%Tbl = Tbl1; 
Tbl = sortrows(Tbl, 'Time');
t = seconds(Tbl.Time); 
xch = 2; x = Tbl{:,xch}; xname = Tbl.Properties.VariableNames{xch}
%x = Tbl.CLFP_NP1___Posterior; xname = 'CLFP_NP1___Posterior';
Fs = 1375; spkFs = 44000; % Hz
load('/Users/torenarginteanu/Desktop/Data_PD/PD24N009/Neuro Omega/mat/LT/Saved To Table/Table Data 2026-09-04 16.02.54 d0a6723833db66/SpkLt_p1560_1A.mat')
spkTbl = Tbl1A; spkTbl.Properties.VariableNames{xch}
xx = spkTbl{:,xch};
%xx = spkTbl.CSPK_NP1___Posterior; 
tsel = (t >= seconds(spkTbl.Time(1))) & (t <= seconds(spkTbl.Time(end)));
%tsel = (t >= 5700) & (t <= 5800);
%tsel = (t >= 8059) & (t <= 8091);
%tsel = (t >= 8077) & (t <= 8089);
%tsel = (t >= 8064) & (t <= 8074);
%tsel = (t >= seconds(spkTbl.Time(1))) & (t <= 8450);
x = x(tsel); t = t(tsel);

%% time-freq; look for beta bursts 

% unprocessed 
figure; spectrogram(x,1*Fs,[],[],Fs,"yaxis","power"); ylim([0 200]);
title(xname);
figure; cwt(x,"bump",Fs);

% adjusted 
[S,fS,tS] = spectrogram(x,1*Fs,[],[],Fs,"yaxis","power");
[~,k1,c2] = pinkcorrect(mean(abs(S),2),fS);
Anoise = k1*fS.^c2; Anoise(1)=eps;
SS = abs(S)./Anoise;
figure; img = imagesc(tS, fS(2:end), (SS(2:end,:))); %colorbar
img.Parent.YDir = 'normal';
title([xname,' adjusted spectrogram']);
ylabel('Frequency (Hz)'); xlabel('time (s)');

%% load spike-sorted data 
load('/Users/torenarginteanu/Desktop/Data_PD/PD24N009/Neuro Omega/mat/LT/Saved To Table/Table Data 2026-09-04 16.02.54 d0a6723833db66/crudeSpikeSort3/times_waveclusdata_LT1A_2.mat')
tSpk = cluster_class(:,2)/1000 + seconds(spkTbl.Time(1));
kidx = cluster_class(:,1);
ku = unique(kidx);
ku = ku(ku > 0); % exclude noise
tSpkK = cell(size(ku));
for ki = 1:length(ku)
    tSpkK{ki} = tSpk(kidx == ku(ki));
    tSpkK{ki} = ...
        tSpkK{ki}(tSpkK{ki}>=t(1) & tSpkK{ki}<=t(end)); % limit spk/lfp to same time selection
end
tSpk = tSpk(kidx > 0); % exclude noise
tSpk = tSpk(tSpk>=t(1) & tSpk<=t(end)); % limit spk/lfp to same time selection

% show spike detection 
figure; plot(spkTbl.Time, xx); hold on; grid on;
for ki = 1:length(ku)
    iSpkKi = round( cluster_class(kidx==ki,2) * spkFs/1000 );
    plot(spkTbl.Time(iSpkKi), xx(iSpkKi), 'o');
end
title('Spike signal detection');
legend(["Signal", "cluster " + (1:length(ku))]);

%% raster plot 

figure; ax(1) = subplot(2,1,1);
if length(ku) > 1
    plot(tSpk,zeros(size(tSpk)),'|', 'LineWidth',2, 'MarkerSize',10); hold on;
    for ki = 1:length(ku)
        tSpkKi = tSpkK{ki};
        plot(tSpkK{ki}, ki*ones(size(tSpkK{ki})), '|', 'LineWidth',2, 'MarkerSize',10);
    end
    ylabel('Spike raster');
end

ax(2) = subplot(2,1,2);
hpf = fir1(2048,minFreq/(Fs/2),"high");
xh = filtfilt(hpf,1,x); xl = x-xh;
plot(t, xh); grid on; 
xlabel('time (s)'); ylabel('LFP');
linkaxes(ax,'x');

%% build spike train and get ISI

[z,dt,Fu,Tl,OnRate,OnOffThresh] = getISIperiods(tSpk, Fs, t);
disp("All Spikes");
disp(" - Total rate: "+num2str(length(tSpk)/(t(end)-t(1))))
disp(" - Mean rate: "+num2str(mean(z)))
disp(" - On Rate: "+num2str(OnRate))

figure; hh = histogram(dt, 0:.001:0.1); grid on; hold on;
hs = stem([Tl, 1/Fu], max(hh.Values)*[1,1], '.');
hs.SeriesIndex = hh.SeriesIndex;
title('ISI histograms:')
subtitle('all spikes + each cluster');
xlabel('ISI (s)'); ylabel('count');
% perform smoothing 
zw  = smoothdata(z,1,'gaussian', ceil(Fs*Tl));
zww = smoothdata(z,1,'gaussian', ceil(Fs*2/Fu));
zw2 = smoothdata(z,1,'gaussian', ceil(1*Fs));

zk = cell(size(ku)); Fuk = nan(size(ku)); Tlk = nan(size(ku));
OnRates = nan(size(ku)); OnOffThreshs = nan(size(ku));
zkw = cell(length(ku),2);
for ki = 1:length(ku)
    tSpkKi = tSpkK{ki};
    [zk{ki},dtki,Fuk(ki),Tlk(ki),OnRates(ki),OnOffThreshs(ki)] = ...
        getISIperiods(tSpkKi, Fs, t);

    disp("Cluster "+num2str(ki))
    disp(" - Total rate: "+num2str(length(tSpkKi)/(t(end)-t(1))))
    disp(" - Mean rate: "+num2str(mean(zk{ki})))
    disp(" - On Rate: "+num2str(OnRates(ki)))

    hh = histogram(dtki, 0:.001:0.1); 
    hs = stem([Tlk(ki), 1/Fuk(ki)], max(hh.Values)*[1,1], '.');
    hs.SeriesIndex = ki+1;
    hh.SeriesIndex = ki+1;
    % perform smoothing
    zkw{ki,1} = smoothdata(zk{ki},1,'gaussian', ceil(Fs*Tlk(ki)));
    zkw{ki,2} = smoothdata(zk{ki},1,'gaussian', ceil(Fs*2/Fuk(ki)));
    zkw{ki,3} = smoothdata(zk{ki},1,'gaussian', ceil(Fs*1));
end

%% report mod components 

TimeFreqPlot(Fs, x, z,zw,zw2, Tl,Fu,OnOffThresh,OnRate, [minFreq, Inf], SS,fS,tS);
sgtitle('All Spikes');

for ki = 1:length(ku)
    try
    TimeFreqPlot(Fs,x, zk{ki},zkw{ki,1},zkw{ki,3}, ...
        Tlk(ki),Fuk(ki),OnOffThreshs(ki),OnRates(ki), ...
        [minFreq, Inf], SS,fS,tS);
    sgtitle(['Spike Cluster ',num2str(ki)]);
    catch ME
        warning(['Issue with cluster ',num2str(ki),': ',ME.message]);
    end
end

%% report rolling rate offset, phase, and amp

Tf = 2*ceil(Tl);
zf = zerocrossrate(zw-zw2, 'method','difference', ...
    'WindowLength',Tf, 'OverlapLength',Tf-1)*Fs/2;
tf = t((Tf/2):(end-Tf/2));
amp2 = envelope(zw-zw2);
figure; sgtitle('All Spikes')
ax2(1) = subplot(2,1,1); 
patch([t; flipud(t)], [zw2; flipud(zw2)]+[amp2; -flipud(amp2)], 'b', ...
    'FaceAlpha',0.5, 'EdgeColor','none');
hold on; grid on; plot(t, zw2, 'b', 'LineWidth', 2);
ylabel('rate (Hz)'); xlabel('t (s)');
title('offset and amplitude');
ax2(2) = subplot(2,1,2);
plot(tf, zf); grid on; 
ylabel('rate (Hz)'); xlabel('t (s)');
title('modulation rate');
linkaxes(ax2, 'x');
xlim([t(1), t(end)]);

for ki = 1:length(zk)
Tf = 2*ceil(Tlk(ki));
zf = zerocrossrate(zkw{ki,2}-zkw{ki,3}, 'method','difference', ...
    'WindowLength',Tf, 'OverlapLength',Tf-1)*Fs/2;
tf = t((Tf/2):(end-Tf/2));
amp2 = envelope(zkw{ki,2} - zkw{ki,3});
figure; sgtitle(['Spike Cluster ',num2str(ki)])
ax2(1) = subplot(2,1,1); 
patch([t; flipud(t)], [zkw{ki,3}; flipud(zkw{ki,3})]+[amp2; -flipud(amp2)], 'b', ...
    'FaceAlpha',0.5, 'EdgeColor','none');
hold on; grid on; plot(t, zkw{ki,3}, 'b', 'LineWidth', 2);
ylabel('rate (Hz)'); xlabel('t (s)');
title('offset and amplitude');
ax2(2) = subplot(2,1,2);
plot(tf, zf); grid on; 
ylabel('rate (Hz)'); xlabel('t (s)');
title('modulation rate');
linkaxes(ax2, 'x');
xlim([t(1), t(end)]);
end

%% compare spike and LFP signals 

% normalize 
xn = (x-mean(x))/std(x);
zkwn = zkw(:,1); zwn = zw; 
%zkwn = cellfun(@(zi) (zi-mean(zi))/std(zi), zkwn, 'UniformOutput',false);
%zwn = (zwn-mean(zwn))/std(zwn);

% evaluate 
r = corr(xn, zwn); 
rr = cellfun(@(zi) corr(xn,zi), zkwn);
lgd = string(ku)+": \rho="+string(rr);
lgd = ["raw LFP"; lgd; "all spk: \rho="+string(r)];

% time domain 
figure; 
plot(t, xn*std(zw) + mean(zwn)); hold on; grid on;  
if length(ku) > 1
for ki = 1:length(zkwn)
    plot(t, zkwn{ki}); % Plot each spike train for comparison
end
end
plot(t, zwn, 'w');
title('Spike-LFP time domain comparison');
xlabel('time (s)'); ylabel('normalized LFP/count');
legend(lgd);

%% see Time-Freq of both LFP and SPK

figure; 
ax3(1) = subplot(2,1,1); 
[W,fW] = cwt(x,"amor",Fs); W = abs(W); fW = log10(fW);
imagesc(t,fW,W); colorbar;
ax3(1).YDir = 'normal';
title('LFP scalogram');
ylabel('log_{10}( freq [Hz] )');
ax3(2) = subplot(2,1,2); 
[W,fW] = cwt(zw,"amor",Fs); W = abs(W); fW = log10(fW);
imagesc(t,fW,W); colorbar;
ax3(2).YDir = 'normal';
title('All Spikes Scalogram');
ylabel('log_{10}( f [Hz] )');
linkaxes(ax3, 'xy');
ylim([-1, 2.3]);
xlabel('time (s)')

for ki = 1:length(zkwn)
figure; 
ax3(1) = subplot(2,1,1); 
[W,fW] = cwt(x,"amor",Fs); W = abs(W); fW = log10(fW);
imagesc(t,fW,W); colorbar;
ax3(1).YDir = 'normal';
title('LFP scalogram');
ylabel('log_{10}( freq [Hz] )');
ax3(2) = subplot(2,1,2); 
[W,fW] = cwt(zkw{ki,1},"amor",Fs); W = abs(W); fW = log10(fW);
imagesc(t,fW,W); colorbar;
ax3(2).YDir = 'normal';
title(['Cluster ',num2str(ki),' Spk Scalogram']);
ylabel('log_{10}( f [Hz] )');
linkaxes(ax3, 'xy');
ylim([-1, 2.3]);
xlabel('time (s)')
end

%% time domain v2 
figure; 
plot(t, xn*(mean(1./diff(tSpk)))); hold on; grid on;  
if length(ku) > 1
for ki = 1:length(ku)
    tSpkKi = tSpkK{ki};
    plot(tSpkKi(2:end), 1./diff(tSpkKi)); % Plot each spike train for comparison
end
end
plot(tSpk(2:end), 1./diff(tSpk), 'w');
title('Spike-LFP time domain comparison');
xlabel('time (s)'); ylabel('normalized LFP/rate');
lgd = ["raw LFP"; string(ku); "all spk"];
legend(lgd);

%% freq domain 
[px,f] = pwelch(xn,[],[],[],Fs,'power'); 
[pz,f] = pwelch(z,[],[],[],Fs,'power');
px = pinkcorrect(px,f); pz = pinkcorrect(pz,f);
px = 20*log10(px); pz = 20*log10(pz);
figure; plot(f,px); hold on; grid on; 
xlabel('freq (Hz)'); ylabel('Power (dB)');
if length(ku) > 1
for ki = 1:length(zk)
    [pzi,f] = pwelch(zk{ki},[],[],[],Fs,'power');
    %pzi = pinkcorrect(pzi,f); pzi = 20*log10(pzi);
    pzi = 20*log10(pzi);
    plot(f,pzi);
end
end
plot(f,pz, 'w'); 
title('(adjusted) spectrum')

%% modulated pulse train analysis 
%{

iwin = (1:(5*Fs)) + 10*Fs;
xwin = xn(iwin); zwin = z(iwin); twin = t(iwin);

[px,frange] = pwelch(xwin,[],[],[],Fs,'power');
X = complex(zeros(size(frange)));
addCompAtFreq = @(fnew,Y,fY,fX) interp1(fY+fnew, Y, fX);

fc = length(spkIdx)/(t(end)-t(1)); % initial est carrier freq = avg spk rate 
zw = movsum(zwin,Fs); % est inst freq
zw0 = mean(zw); zw = zw-zw0; % different from carrier freq?
r = zw0/fc; % "repetition rate"

% inst freq fourier sine coeffs 
[fsine, ampsine, phsine, amp0] = FourierSine(zw, Fs, twin);

% pulse train fourier coeffs 
wf = mean(WF); %wf = wf-mean(wf); % 0 offset
wfFT = fftshift(fft(wf)); % kernel for conv 
%wfFTres = fs/(length(wfFT)-1); % hz per sample
%wfFT = resample(wfFT,round(Fs),round(fs)); % for conv with z
wff = linspace(-fs/2, fs/2, length(wf));
wfFTresamp = resample(wfFT,round(fs),round(Fs)); % should be same res as frange
wfFTresamp_ihalf = length(wfFTresamp)/2;
wfFTresamp(1:wfFTresamp_ihalf) = conj(fliplr(wfFTresamp(wfFTresamp_ihalf+1:end)));

% k=0 component 
%{
f0 = 2*pi*fsine'; f0 = [-fliplr(f0),0,f0];
q0 = [flipud(ampsine).*exp(-1i*flipud(phsine)); fc; ampsine.*exp(1i*phsine)]';
%}
f0 = [0;  fsine]';
q0 = [fc; ampsine.*exp(1i*phsine)/(2*r)]';
for fi = 1:length(f0)
    %X = X + q0(fi)*addCompAtFreq(f0(fi), wfFT, wff, frange);
    [~,fj] = min(abs(f0(fi)-frange));
    X(fj) = X(fj) + q0(fi);
end

k = 1;
while k*fc < 120
    [fk,qk] = getFQ(k, r, fc, fsine, ampsine, phsine);
    for fi = 1:length(fk)
        %X = X + qk(fi)*addCompAtFreq(fk(fi), wfFT, wff, frange);
        [~,fj] = min(abs(fk(fi)-frange));
        X(fj) = X(fj) + qk(fi);
    end
    k = k+1;
end

figure; 
plot(frange, 20*log10(px)); hold on; grid on; 
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
plot(frange, 20*log10(abs(X).^2), '.');

%}
%% helpers

function [z,dt,Fu,Tl,OnRate,OnOffThresh]...
    = getISIperiods(tSpk, Fs, t)

% build pulse train
t0 = t(1);
dthalf = 1/(2*Fs);
z = arrayfun(@(ti) sum((tSpk>(ti-dthalf))&(tSpk<(ti+dthalf))), t)/(2*dthalf);
iSpk = (tSpk-t0)*Fs; iSpk = round(iSpk) + 1;

% average ("carrier") firing rate in 1s sliding windows: 
zw2 = smoothdata(z,1,'gaussian', ceil(1*Fs));

% distinguish on/off states 
[~,OnOffThresh] = midcross(zw2); 
OnOffRate = statelevels(zw2); OnRate = OnOffRate(2);
isOff = zw2 < OnOffThresh;
isOff = movmean( double(isOff), 5 ) > 0;
spkOn = ~isOff(iSpk);

% ISI periods 
dt = diff(tSpk);
dt = dt(spkOn(2:end) & spkOn(1:(end-1)));
%[~,~,Fu] = isoutlier(1./dt, "ThresholdFactor",3); % upper bound on carrier freq
%[~,~,Tl] = isoutlier(dt,    "ThresholdFactor",1); % 1/ lower bound on carrier freq
%[~,Tu,Tl] = isoutlier(dt, "quartiles", "ThresholdFactor",0.5); 
[~,Tu,Tl] = isoutlier(dt, "percentiles", [1 80]); 
Fu = 1/Tu;
end


function [F,W,fW,Wt] = getModComponents(zw2, Fs, OnOffThresh, fminmax)
if nargin < 4
    fminmax = [];
    if nargin < 3
        OnOffThresh = 0;
    end
end
if isempty(fminmax)
    fminmax = [0 inf];
end

% cont wavelet transform pks = candidate components 
[W,fW] = cwt(zw2,"amor",Fs);
iOn = zw2 > OnOffThresh; W = W(:, iOn);
Wt = mean(abs(W),2);
[vPk,fPk,wPk,pPk] = findpeaks(flipud(Wt), flipud(fW), "WidthReference","halfheight");

% exclude freqs too low/high 
selPk = (fPk >= fminmax(1)) & (fPk <= fminmax(2));
vPk = vPk(selPk); 
fPk = fPk(selPk); 
wPk = wPk(selPk); 
pPk = pPk(selPk); 

% exclude small pks 
selPk = (pPk >= max(Wt)/10) | ((pPk >= max(Wt)/200) & (vPk >= max(Wt)*.4));
vPk = vPk(selPk); 
fPk = fPk(selPk); 
wPk = wPk(selPk); 
pPk = pPk(selPk); 

% adjust/report out freq ranges 
F = [fminmax(1).*[1 1 1]; 
     fPk + [-.5,0,.5].*wPk;  
     fminmax(2).*[1 1 1] ];
for iPk = 2:height(F)-1
    F(iPk,1) = max(F(iPk,1), F(iPk-1,3));
end
for iPk = fliplr(2:height(F)-1)
    F(iPk,3) = min(F(iPk,3), F(iPk+1,1));
end
F = F(2:(end-1),:);

end


function [A, k1, c2] = pinkcorrect(A,f)
% correct for noise that obeys Anoise = k1*f^c2
% i.e. ln(Anoise) = c2*ln(f) + c2*ln(k1)
if f(1) < 2*eps
    f0 = 0; f = f(2:end);
    A0 = A(1,:); A = A(2:end,:);
else
    f0 = zeros(0,width(f)); A0 = zeros(0,width(A));
end
lnA = log(A); lnf = log(f); F = [ones(size(lnf)), lnf];
c = F\lnA; 
% c1 = c2*ln(k1), i.e. k1 = exp(c1/c2)
c2 = c(2); k1 = exp(c(1)/c(2));
lnAnoise = F*c;
lnA = lnA - lnAnoise; A = exp(lnA);
A = [A0; A];
end


function [fk,qk] = getFQ(k, r, fc, fsine, ampsine, phsine)
% for k > 0 only! 
% positive sided only!
maxnumel = 1e9;

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
maxnumeladj = (maxnumel^(1/length(ampsine)))/(2); 
maxnumeladj = floor(maxnumeladj);
Jmag = sort(J(:).^2, 'descend'); Jmin = Jmag(maxnumeladj);
powerloss = sum(Jmag((maxnumeladj+1):end))/sum(Jmag);
warning(['Losing ',num2str(100*powerloss,2),'% power.'])
J(J.^2 < Jmin) = nan;

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
    if ~isempty(ql)
        try
        ff1 = (fl'+f1); f1 = single(ff1(:))';
        qq1 = ql'*q1; q1 = single(qq1(:))';
        catch ME
            keyboard
        end
    end
end
f1=f1+k*fc; q1 = q1*fc;

% second term (freq shifted)
fshift = [-flipud(fsine);fsine]';
Q1b = [flipud(ampsine).*exp(-1i*flipud(phsine)); ampsine.*exp(1i*phsine)]';
Q1b = reshape(Q1b,[1,1,length(Q1b)]); fshift = reshape(fshift,[1,1,length(fshift)]);
Q1b = Q1b.*Q1; F1b = (fsine*nrange)' + fshift;
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
        try
        ff1 = (fl+f1b); f1b = single(ff1(:))';
        qq1 = ql'*q1b; q1b = single(qq1(:))';
        catch MEb
            keyboard
        end
    end
end
f1b=f1b+k*fc; q1b = q1b/(2*r);

% output 
fk = [f1, f1b]; qk = [q1, q1b];

end


function [fsine, ampsine, phsine, A0] = FourierSine(zw, Fs, t, K, tres, fres, fminmax)
if nargin < 7
    fminmax = [];
    if nargin < 6
        fres = false;
        if nargin < 5
            tres = false;
            if nargin < 4
                K = 8;
                if nargin < 3
                    t = ((1:length(zw))-1)/Fs;
                end
            end
        end
    end
end

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

if fres
    fresnorm = ceil(fres*Nt/Fs); % freq samples
end

% build sine-wave-only representation: A_k * sin(2*pi*f_k*t + phi_k)
if tres 
    [S,fS,tS] = spectrogram(y,tres*Fs,[],[],Fs,"yaxis");
    S = S/(2*Fs);
    Ak = mean(abs(S), 2);
    phk = circ_mean(angle(S), [], 2);
    fk = fS;
else
    Ak = Apos;
    fk = freqs(posIdx);
    phk = phpos - pi/2; % convert complex exponential phase to sine-phase
end
if fres
    Ak = smoothdata(Ak,1,"gaussian", fresnorm);
end

% exclude low freqs 
if ~isempty(fminmax)
    fsel = (fk > fminmax(1)) & (fk < fminmax(2));
    Ak = Ak(fsel);
    phk = phk(fsel);
    fk = fk(fsel);
end

% (since e^{i(ωt+ph)} = cos(ωt+ph) + i sin(...); to express as sin(ωt+phi_s):
% sin(ωt+phi_s) = cos(ωt+phi_s - pi/2), adjust accordingly. Using phase shift:
% amplitude already accounts for factor 2 above.)

% reconstruct using a limited number of harmonics for stability (optional)
% choose K components with largest amplitudes
K = min(K, numel(Ak)); % limit to _ components or less
if fres
    [pkA,pkF,pkFw,pkAp] = findpeaks(Ak, 'SortStr','descend', ...
        'MinPeakDistance',fresnorm); % could instead use 'MinPeakWidth' and sort by Ap
    K = min(K, numel(pkF));
    sel = pkF(1:K);
else
    [~,ord] = sort(Ak,'descend');
    sel = sort(ord(1:K));
end

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


function TimeFreqPlot(Fs, x, z,zw,zw2, Tl,Fu,OnOffThresh,OnRate, MinMaxModFreq, SS,fS,tS)

t = 1:length(x); t = (t-1)/Fs; t = t';

% Multi-panel spectrogram with linked marginal plots
[F,W,fW,Wmean] = getModComponents(zw, Fs, OnOffThresh, MinMaxModFreq);
fPlot = fS(2:end);
sPlot = SS(2:end,:);
specMean = mean(sPlot, 2, 'omitnan');
BPFc = fir1(1023, [1/Tl, Fu/2]/(Fs/2)); 
%BPFc = fir1(1023, (OnRate+[-5,5])/(Fs/2)); 
xc = filtfilt(BPFc,1,x); zc = filtfilt(BPFc,1,z); % carrier band

% manually enforce a max component band width 
F(:,3) = min(F(:,3), F(:,2)+5);
F(:,1) = max(F(:,1), F(:,2)-5);

figure;

% main (top-right): time+frequency of LFP 
axMain = axes('Position',[0.25 0.50 0.70 0.45]);
imagesc(axMain, tS, fPlot, sPlot);
axMain.YLim = [0 200];
axMain.YDir = 'normal';
ylabel(axMain,'Frequency (Hz)');
title(axMain,'LFP Time-Freq. Analysis (adjusted spectrogram)');
colormap(axMain, parula);
%colorbar(axMain);

% left: frequency analysis of spike train 
axLeft = axes('Position',[0.06 0.50 0.14 0.45]);
h1 = histogram(axLeft, zw(zw>MinMaxModFreq(1)), 'Orientation','horizontal');
hold(axLeft,'on'); grid(axLeft,'on');
h2 = histogram(axLeft, zw2(zw2>OnOffThresh), 'Orientation','horizontal');
%h2 = histogram(axLeft, zw(zw>MinMaxModFreq(1)),   'Orientation','horizontal');
axLeft_xmax = max(h1.Values(2:end));
axLeft_xmax = max(axLeft_xmax, max(h2.Values(2:end)));
Wmean = Wmean * axLeft_xmax / max(Wmean); % scale
plot(axLeft,[0,axLeft_xmax],OnRate*[1,1], ...
    'LineWidth',3);
plot(axLeft,Wmean,fW, ...
    'LineWidth',2); 
%{
errorbar(axLeft, axLeft_xmax*ones(height(F),1), ...
    F(:,2), F(:,1)-F(:,2), F(:,3)-F(:,2), [],[], ...
    'LineWidth', 2);
%}
legend(axLeft, ...
    ['Inst Fire Rate ',num2str(Tl),'s'], ...
    'Avg Fire Rate 1s', ...
    'Est ''On'' Fire Rate', ...
    'Fire Rate Mod. Spectrum', ...
    ...'Sel. Components', ...
    'Location','best');
axLeft.XDir = "reverse";
axLeft.YDir = 'normal';
axLeft.YLim = axMain.YLim;
axLeft.XLim = [0, 1.1*axLeft_xmax];
axLeft.YAxisLocation = 'right';
xlabel(axLeft,'Power / count');
title(axLeft, 'Spike Rate Freq. Analysis');
axLeft.XTick = [];
linkaxes([axLeft axMain],'y');

axBottom1 = axes('Position',[0.25 0.33 0.70 0.13]);
plot(axBottom1, t, (xc-mean(xc))/std(xc)); 
grid(axBottom1,'on'); hold(axBottom1,'on');
plot(axBottom1, t, (x-mean(x))/std(x)); 
for Fi = fliplr(1:height(F))
    BPF = buildFIRBPF(Fs,F(Fi,1),F(Fi,3),2,1023);
    xFi = filtfilt(BPF,1,x);
    plot(axBottom1, t, (xFi-mean(xFi))/std(xFi));
end
ylabel(axBottom1, 'LFP (norm.)');
legend(axBottom1, ["carrier"; "full"; string(flipud(F(:,2)))+" Hz"]);
title(axBottom1, 'Time Domain Analysis');
axBottom1.XLim = axMain.XLim; axBottom1.XTickLabel = {};

axBottom2 = axes('Position',[0.25 0.19 0.70 0.13]);
amp2 = envelope(zw-zw2);
plot(axBottom2, t, zw2, 'b', 'LineWidth', 2);
grid(axBottom2,'on'); hold(axBottom2,'on');
patch(axBottom2, [t; flipud(t)], [zw; flipud(zw)]+[amp2; -flipud(amp2)], 'b', ...
    'FaceAlpha',0.5, 'EdgeColor','none');
ylabel(axBottom2, 'Fire Rate');
legend(axBottom2, 'Avg', 'Mod Amp');
axBottom2.XLim = axMain.XLim; axBottom2.XTickLabel = {};

axBottom3 = axes('Position',[0.25 0.05 0.70 0.13]);
grid(axBottom3,'on'); hold(axBottom3,'on');
for Fi = 1:height(F)
    BPF = buildFIRBPF(Fs,F(Fi,1),F(Fi,3),2,1023);
    xFi = filtfilt(BPF,1,x); zFi = filtfilt(BPF,1,z);
    [xPAC,txPAC] = calcPACwindowed(xFi,xc,10*Fs,10);
    [zPAC,tzPAC] = calcPACwindowed(zFi,zc,10*Fs,10);
    hx = plot(t(txPAC), xPAC, 'LineWidth',1);
    hx.SeriesIndex = Fi;
    hz = plot(t(tzPAC), zPAC, '--', 'LineWidth',1);
    hz.SeriesIndex = Fi;
end
lgd = F(:,2)+" Hz "+["LFP", "SPK"]; lgd = lgd'; lgd = lgd(:);
legend(axBottom3, lgd);
ylabel(axBottom3, 'PAC (MI)');
xlabel(axBottom3, 'Time (s)');
axBottom3.XLim = axMain.XLim; 
linkaxes([axMain axBottom1 axBottom2 axBottom3],'x');

end