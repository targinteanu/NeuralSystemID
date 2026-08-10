%% PY26N010 analysis
% phase-dependent stim and CCEPs

%% phase-dep CCEPs 

% data of interest 
nstbl = ns2timetable('/Users/torenarginteanu/Desktop/Data_EMU/PY26N010/Phase-Dependent CCEPs/PDCCEP001.ns2');
X = [nstbl.RA1, nstbl.RAH2]; 

% timing 
fs = nstbl.Properties.SampleRate;
t = nstbl.Time; 
t0 = t(1); 
T = [minutes(8),    minutes(9);    ... pk 66-65-67 (RAH2-1-3) 3mA 0.5hz  
     minutes(9),    minutes(10);   ... ibid thresh100
     minutes(14),   minutes(16.5); ... pk 1(RA1)-65-66 4mA 0.5hz 40+trl thresh100
     minutes(16.5), minutes(19)];    % tr 66-65-67 4mA 0.5hz 40+trl thresh100
t1 = datetime(2026,8,5,19,0,0); t1.TimeZone = t0.TimeZone; T = T+t1;
tRel = seconds(t-t0);

% stim timing 
g = nstbl.AINP1; g = [false; diff(g) > 1000];
gidx = find(g); 

% art remove 
artdur = 400; % samples 
G = movmean(g,artdur)>eps;
Gshift = floor(artdur/2) - 1;
G = [false(Gshift,1); G(1:(end-Gshift))];
tRelNoArt = tRel(~G);
XNoArt = X(~G,:);
XArtRem = interp1(tRelNoArt, XNoArt, tRel, 'linear', 'extrap');

% filter, phase 
BPF = buildFIRBPF(fs,4,9);
Xf = filtfilt(BPF,1,XArtRem);
Xph = angle(hilbert(Xf));
figure; 
for p = 1:width(X)
    ax(p) = subplot(width(X),1,p);
    plot(t, X(:,p)-mean(XArtRem(:,p))); hold on; grid on; 
    plot(t, XArtRem(:,p)-mean(XArtRem(:,p)));
    plot(t, Xf(:,p))
end
linkaxes(ax, 'x');

% analysis 
figure; getPhaseEP(X(:,2),Xph(:,2),g,t,[T(1,1),T(2,2)],[-250 749],fs);
sgtitle('RAH2 peak 3mA stim RAH1-RAH3'); 
figure; getPhaseEP(X(:,2),Xph(:,2),g,t,[T(4,1),T(4,2)],[-250 749],fs);
sgtitle('RAH2 trough 4mA stim RAH1-RAH3'); 
figure; getPhaseEP(X(:,1),Xph(:,1),g,t,[T(3,1),T(3,2)],[-250 749],fs);
sgtitle('RA1 peak 4mA stim RAH1-RAH2'); 

%% phase-target stim 
% rec channel = RAH2 
% encode: peak target; decode: trough target 

nstbl = ns2timetable('/Users/torenarginteanu/Desktop/Data_EMU/PY26N010/Memory/memory001.ns2');
evtbl = parseNEVSerialCodes('/Users/torenarginteanu/Desktop/Data_EMU/PY26N010/Memory/memory001.nev');

tstart = evtbl(strcmp(evtbl.EventName, 'decoding_start'),:);
tend = evtbl(strcmp(evtbl.EventName, 'decoding_end'),:);
tend = tend(1:2:end,:); % remove duplicates 
Tdecode = [tstart.Time, tend.Time]; 

tstart = evtbl(strcmp(evtbl.EventName, 'image_1_adaptive'),:);
tend = evtbl(strcmp(evtbl.EventName, 'end_encoding'),:);
Tencode = [tstart.Time, tend.Time]; 

%% helper(s)
function getPhaseEP(x,xph,trig,t,trng,EPidx,fs)
tEP = EPidx(1):EPidx(2); tEP = tEP/fs;
xidx = (t > trng(1)) & (t < trng(2));
x = x(xidx); xph = xph(xidx); trig = find(trig(xidx));
EPl = diff(EPidx)+1;
EP = nan(EPl,length(trig));
for ti = 1:length(trig)
    EP(:,ti) = x(trig(ti)+(EPidx(1):EPidx(2)));
end
bldur = 5-EPidx(1);
EPos = EP(1:bldur,:); EPos = mean(EPos);
EP_ = EP - EPos; 
EPm = mean(EP_,2); EPv = std(EP_,[],2);
ph = xph(trig); phmean = circ_mean(ph); phconf = circ_std(ph);
phmean = phmean*180/pi; phconf = phconf*180/pi;
subplot(2,1,1); polarhistogram(ph, 12);
title('Stim Phase');
subtitle([num2str(phmean),'±',num2str(phconf),'°']);
subplot(2,1,2); 
errorbar(tEP,EPm,EPv); grid on;
xlabel('sec'); ylabel('uV'); title('Evoked Potential');
end