ViewSpectrum
%%
t = dtaBL.Time; x = dtaBL.(channelNameRec);
t = seconds(t-t(1));

BPFbeta = buildFIRBPF(SampleRate, 13,30, 2, 201);
BPFglo = buildFIRBPF(SampleRate, 30,70, 2, 201);
BPFghi = buildFIRBPF(SampleRate, 70,150, 2, 201);
qFactor = 35;
[n60b, n60a] = iirnotch(60/(SampleRate/2), (60/(SampleRate/2))/qFactor);
[n120b, n120a] = iirnotch(120/(SampleRate/2), (120/(SampleRate/2))/qFactor);

xB = filtfilt(BPFbeta,1,x); 
xGlo = filtfilt(n60b,n60a,x); xGlo = filtfilt(BPFglo,1,xGlo);
xGhi = filtfilt(n120b,n120a,x); xGhi = filtfilt(BPFghi,1,xGhi); 

%% time-frequency 

[S,fS,tS] = spectrogram(x, [],[],[], SampleRate);
S = abs(S); % magnitude 
S = 20*log10(S); % power (dB)
S = flipud(S); % low f at bottom 
%[p,f] = periodogram(x,[],[],SampleRate,'power');
[p,f] = pwelch(x,[],[],[],SampleRate,'power');
p = 10*log10(p); % power (dB)

figure; sgtitle(channelNameRec);

ax(1,1) = subplot(2,2,1);
imagesc(S); colorbar; 
ylabel('freq');
yticklabels(fS(fliplr(yticks))); xticklabels(tS(xticks));

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

%linkaxes(ax(1,:), 'y'); 
%linkaxes(ax(:,1), 'x');

%% amp-phase 