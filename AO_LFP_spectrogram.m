%% load raw data 

[fn,fp] = uigetfile('*.mat');
load(fullfile(fp,fn));
Tbl = sortrows(Tbl, 'Time');
t = seconds(Tbl.Time);
Fs = 1/median((diff(t))); % Hz
tReg = t(1):(1/Fs):t(end);

E = Tbl.Properties.Events;
E = E(contains(E.EventLabels, 'Stim'), :);
D = AlphaOmegaTable2Depth(Tbl);

%{
% truncate time to match spkTbl
tsel = (t >= (spkTbl.Time(1))) & (t <= (spkTbl.Time(end)));
t = t(tsel);
Tbl = Tbl(tsel,:); 
Esel = (E.Time >= (spkTbl.Time(1))) & (E.Time <= (spkTbl.Time(end)));
E = E(Esel,:);
Dsel = (D.Time >= (spkTbl.Time(1))) & (D.Time <= (spkTbl.Time(end)));
D = D(Dsel,:);
%}

%%
for xch = 1:width(Tbl)
    try

        x = Tbl{:,xch}; xname = Tbl.Properties.VariableNames{xch};
        xname(xname=='_') = ' ';

% regularize 
[t,iu] = unique(t); x = x(iu);
x = interp1(t,x,tReg, "nearest","extrap");

% plug NaNs 
inan = isnan(x); 
if any(inan)
    tNan = tReg(~inan); xNan = x(~inan);
    x = interp1(tNan,xNan,tReg, "nearest","extrap");
end

%% spectrogram; look for beta bursts 

% compute spectrogram
%figure; spectrogram(x,1*Fs,[],[],Fs,"yaxis","power"); ylim([0 200]);
%title(xname);
[S,fS,tS] = spectrogram(x,1*Fs,[],[],Fs,"yaxis","power");
tS = tS + tReg(1);

% correct pink noise
[~,k1,c2] = pinkcorrect(mean(abs(S),2),fS);
Anoise = k1*fS.^c2; Anoise(1)=eps;
SS = abs(S)./Anoise;

% saturate out outliers for better display
%SSall = log(SS(:)+eps);
SSall = SS(:);
[~,~,OLthresh] = isoutlier(SSall, 'median', 'ThresholdFactor',10);
OLthresh = max(SSall(SSall<OLthresh));
%OLthresh = exp(OLthresh);

% plot spectrogram, LFP, depth, stim
figure; 
ax(1) = subplot(311);
img = imagesc(tS, fS(2:end), (SS(2:end,:)), [0,OLthresh]); %colorbar
img.Parent.YDir = 'normal';
title([xname,' adjusted spectrogram']);
ylabel('Frequency (Hz)'); xlabel('time (s)');
ax(2) = subplot(312);
plot(tReg, x); grid on; ylabel('LFP sig')
ax(3) = subplot(313);
if ~isempty(D)
plot(seconds(D.Time), D.DEPTH); ylabel('Depth (mm)');
end
hold on; grid on; 
if ~isempty(E)
stem(seconds(E.Time), -5*ones(size(E)), '.'); 
end
if ~isempty(D) && ~isempty(E)
    legend('Depth', 'Stim');
end
linkaxes(ax, 'x');
    catch ME
        warning(ME.message)
    end
end