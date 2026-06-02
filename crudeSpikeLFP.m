%% load data 
load('/Users/torenarginteanu/Desktop/Data_PD/PD26N003/Neuro Omega/SavedTable1375HzRT.mat')
Tbl = sortrows(Tbl, 'Time');
t = seconds(Tbl.Time);
x = Tbl.CLFP_AP_T___Central;
Fs = 1375; % Hz
load("/Users/torenarginteanu/Desktop/Data_PD/PD26N003/Neuro Omega/SPK_RT_SelectedTimes.mat")
load('/Users/torenarginteanu/Desktop/Data_PD/PD26N003/Neuro Omega/spikesort_227669d9c8fe69f115266d7020ddc4244c5f3779.mat')
spkTbl = depth_p0496_1;
tSpk = spkIdx/fs + seconds(spkTbl.Time(1));
tSpk1 = spkIdx(kidx==1)/fs + seconds(spkTbl.Time(1));
tSpk2 = spkIdx(kidx==2)/fs + seconds(spkTbl.Time(1));
tsel = (t >= seconds(spkTbl.Time(1))) & (t <= seconds(spkTbl.Time(end)));
x = x(tsel); t = t(tsel);

%% filter LFP 
hpf = fir1(1024,0.25/(Fs/2),"high");
xh = filtfilt(hpf,1,x); xl = x-xh;

%% build spike train 
z = zeros(size(t)); z1 = zeros(size(t)); z2 = zeros(size(t));
for zi = tSpk'
    [~,ti] = min(abs(t-zi));
    z(ti) = 1;
end
for zi = tSpk1'
    [~,ti] = min(abs(t-zi));
    z1(ti) = 1;
end
for zi = tSpk2'
    [~,ti] = min(abs(t-zi));
    z2(ti) = 1;
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
z1w = smoothdata(z1,1,'gaussian',w);
z2w = smoothdata(z2,1,'gaussian',w);

%% compare spike and LFP signals 

% normalize 
xn = (x-mean(x))/std(x);
z1wn = (z1w-mean(z1w))/std(z1w);
z2wn = (z2w-mean(z2w))/std(z2w);
zwn = (zw-mean(zw))/std(zw);

% evaluate 
r = corr(xn, zwn); r1 = corr(xn, z1wn); r2 = corr(xn, z2wn);

% time domain 
figure; 
plot(t, xn); hold on; grid on; plot(t, zwn); plot(t, z1wn); plot(t, z2wn);
title('Spike-LFP time domain comparison');
xlabel('time (s)'); ylabel('normalized');
legend('raw LFP', ...
    ['all: \rho=',num2str(r)], ...
    ['1: \rho=',num2str(r1)], ...
    ['2: \rho=',num2str(r2)]);

% freq domain 
figure; 
pwelch(xn,[],[],[],Fs,'power'); hold on;
pwelch(zwn,[],[],[],Fs,'power');
pwelch(z1wn,[],[],[],Fs,'power'); 
pwelch(z2wn,[],[],[],Fs,'power'); 