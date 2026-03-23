ViewSpectrum
%%
channelName = channelNameRec(1);
t = dtaBL.Time; x = dtaBL.(channelName);
t = seconds(t-t(1));

BPFbeta = buildFIRBPF(SampleRate, 20,40, 2, 201);
BPFglo = buildFIRBPF(SampleRate, 40,80, 2, 201);
BPFghi = buildFIRBPF(SampleRate, 80,160, 2, 201);
qFactor = 35;
[n60b, n60a] = iirnotch(60/(SampleRate/2), (60/(SampleRate/2))/qFactor);
[n120b, n120a] = iirnotch(120/(SampleRate/2), (120/(SampleRate/2))/qFactor);

xB = filtfilt(BPFbeta,1,x); 
xGlo = filtfilt(n60b,n60a,x); xGlo = filtfilt(BPFglo,1,xGlo);
xGhi = filtfilt(n120b,n120a,x); xGhi = filtfilt(BPFghi,1,xGhi); 

%% time-frequency 

[S,fS,tS] = spectrogram(x, 1*SampleRate,[],[], SampleRate);
S = abs(S); % magnitude 
S = 20*log10(S); % power (dB)
%[p,f] = periodogram(x,[],[],SampleRate,'power');
[p,f] = pwelch(x,[],[],[],SampleRate,'power');
p = 10*log10(p); % power (dB)

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

nbin = 10;
betabinedge = linspace(20,40,nbin);
globinedge = 2*betabinedge; ghibinedge = 2*globinedge;
betabincent = .5*(betabinedge(2:end) + betabinedge(1:(end-1)));
globincent = .5*(globinedge(2:end) + globinedge(1:(end-1)));
ghibincent = .5*(ghibinedge(2:end) + ghibinedge(1:(end-1)));
betabinind = discretize(fS', betabinedge); 
globinind = discretize(fS', globinedge); 
ghibinind = discretize(fS', ghibinedge); 
Sbeta = S(betabinind,:); Sglo = S(globinind,:); Sghi = S(ghibinind,:);

figure; 
plot(Sbeta(:), Sglo(:), 'v'); hold on; grid on; plot(Sbeta(:), Sghi(:), '^');
xlabel('\beta power'); ylabel('\gamma power'); legend('lo\gamma', 'hi\gamma');
title(channelName+" Harmonic Analysis");

%% amp-phase 