ViewSpectrum
channelName = channelNameRec(1);
%%
t = dtaBL.Time; x = dtaBL.(channelName);
t = seconds(t-t(1));

BPFbeta = buildFIRBPF(SampleRate, 13,30,  2, 201);
BPFglo  = buildFIRBPF(SampleRate, 30,70,  2, 201);
BPFghi  = buildFIRBPF(SampleRate, 70,150, 2, 201);
qFactor = 35;
[n60b, n60a] = iirnotch(60/(SampleRate/2), (60/(SampleRate/2))/qFactor);
[n120b, n120a] = iirnotch(120/(SampleRate/2), (120/(SampleRate/2))/qFactor);

xB = filtfilt(BPFbeta,1,x); 
x = filtfilt(n60b,n60a,x);
xGlo = filtfilt(BPFglo,1,x);
x = filtfilt(n120b,n120a,x); 
xGhi = filtfilt(BPFghi,1,x); 

%% time-frequency 

% spectrum/ogram
[S,fS,tS] = spectrogram(x, 1*SampleRate,[],[], SampleRate);
S = abs(S); % magnitude 
S = 20*log10(S); % power (dB)
[p,f] = periodogram(x,[],[],SampleRate,'psd');
%[p,f] = pwelch(x,[],[],[],SampleRate,'psd');
p = 20*log10(p); % power (dB)

% pink noise correct 
F = [log10(f), ones(size(f))];
pinkcoef = F(2:end,:)\p(2:end,:);
pinkP = F*pinkcoef;
p = p-pinkP; 
FS = [log10(fS), ones(size(fS))];
pinkcoefS = FS(2:end,:)\mean(S(2:end,:),2);
pinkPS = FS*pinkcoefS;
S = S-pinkPS;

figure; sgtitle(channelName);

ax(1,1) = subplot(2,2,1);
imagesc(tS, fS, S); %colorbar; 
ax(1,1).YDir = 'normal'; % low f at bottom
ylabel('freq');

ax(1,2) = subplot(2,2,2);
plot(p,f); grid on;
xlabel('power (dB)');

ax(2,1) = subplot(2,2,3); 
yyaxis right; plot(t, x); ylabel('unfiltered'); hold on; grid on;
yyaxis left;
plot(t, xB); plot(t, xGlo, '--'); plot(t, xGhi, ':');
ylabel('filtered'); 
legend('\beta', 'lo\gamma', 'hi\gamma', 'Raw');
xlabel('time (s)');

linkaxes(ax(1,:), 'y'); 
linkaxes(ax(:,1), 'x');

%% harmonic analysis 

fSbetaInd = (fS >= 13)&(fS < 30);
fSgloInd  = (fS >= 30)&(fS < 70); 
fSghiInd  = (fS >= 70)&(fS < 150); 
[maxSbeta,maxfSbetaInd] = max(S(fSbetaInd,:));
[maxSglo,maxfSgloInd] = max(S(fSgloInd,:));
[maxSghi,maxfSghiInd] = max(S(fSghiInd,:));
maxfSbeta = fS(fSbetaInd); maxfSbeta = maxfSbeta(maxfSbetaInd);
maxfSglo = fS(fSgloInd); maxfSglo = maxfSglo(maxfSgloInd);
maxfSghi = fS(fSghiInd); maxfSghi = maxfSghi(maxfSghiInd);

figure; 
plot(maxfSbeta, maxfSglo, 'v'); hold on; grid on; plot(maxfSbeta, maxfSghi, '^');
xlabel('\beta peak freq (Hz)'); ylabel('\gamma peak freq (Hz)'); 
legend('lo\gamma', 'hi\gamma');
title(channelName+" Harmonic Analysis");

%% amp-phase 
calcPAC(xB, xGlo, 18, true); hold on;
calcPAC(xB, xGhi, 18, gca());