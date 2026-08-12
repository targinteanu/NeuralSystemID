%% PY26N010 analysis
% phase-dependent stim and CCEPs

%% phase-dep CCEPs 

% data of interest 
nstbl = ns2timetable('/Users/torenarginteanu/Desktop/Data_EMU/PY26N010/Phase-Dependent CCEPs/PDCCEP001.ns2');
X = [nstbl.RA1, nstbl.RAH2]; 

%% timing 
fs = nstbl.Properties.SampleRate;
t = nstbl.Time; 
t0 = t(1); 
T = [minutes(8),    minutes(10);    ... pk 66-65-67 (RAH2-1-3) 3mA 0.5hz  
     ...minutes(9),    minutes(10);   ... ibid thresh100
     minutes(14),   minutes(16.5); ... pk 1(RA1)-65-66 4mA 0.5hz 40+trl thresh100
     minutes(16.5), minutes(19)];    % tr 66-65-67 4mA 0.5hz 40+trl thresh100
t1 = datetime(2026,8,5,19,0,0); t1.TimeZone = t0.TimeZone; T = T+t1;
tRel = seconds(t-t0);
g = nstbl.AINP1; g = [false; diff(g) > 1000];
gidx = find(g);

% baseline 
BPF = buildFIRBPF(fs,4,9);
iBL = (t < T(2,1)) & (t > T(1,2));
XBL = X(iBL,:)-mean(X(iBL));
XBL = filtfilt(BPF,1,XBL);
ARord = 100; dsrat = 10;
XBLd = XBL(1:dsrat:end,:);
mdl1 = ar(iddata(XBLd(:,1),[],1/fs),ARord/dsrat,'yw');
mdl2 = ar(iddata(XBLd(:,2),[],1/fs),ARord/dsrat,'yw');
m1 = mdl1.A; m2 = mdl2.A;
mdl1 = zeros(1,ARord+1); mdl1(1:dsrat:end) = m1;
mdl2 = zeros(1,ARord+1); mdl2(1:dsrat:end) = m2;

% remove artifact 
XArtRem = [artremoveAR(mdl1, ARord, X(:,1), gidx, -10, 400), ...
           artremoveAR(mdl2, ARord, X(:,2), gidx, -10, 400)];

% filter, phase 
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
figure; getplotPhaseEP(X(:,2),Xph(:,2),g,t,[T(1,1),T(1,2)],[-250 749],fs);
sgtitle('RAH2 peak 3mA stim RAH1-RAH3'); 
figure; getplotPhaseEP(X(:,2),Xph(:,2),g,t,[T(3,1),T(3,2)],[-250 749],fs);
sgtitle('RAH2 trough 4mA stim RAH1-RAH3'); 
figure; getplotPhaseEP(X(:,1),Xph(:,1),g,t,[T(2,1),T(2,2)],[-250 749],fs);
sgtitle('RA1 peak 4mA stim RAH1-RAH2'); 

%% phase-target stim 
% rec channel = RAH2 
% encode: peak target; decode (recall): trough target 

nstbl = ns2timetable('/Users/torenarginteanu/Desktop/Data_EMU/PY26N010/Memory/memory001.ns2');
evtbl = parseNEVSerialCodes('/Users/torenarginteanu/Desktop/Data_EMU/PY26N010/Memory/memory001.nev');
%%
fs = nstbl.Properties.SampleRate;
t = nstbl.Time; X = nstbl.RAH2;
t0 = t(1); 
tRel = seconds(t-t0);
g = nstbl.AINP1; g = [false; diff(g) > 1000];

tstart = evtbl(strcmp(evtbl.EventName, 'decoding_start'),:);
tend = evtbl(strcmp(evtbl.EventName, 'decoding_end'),:);
tend = tend(1:2:end,:); % remove duplicates 
Tdecode = [tstart.Time, tend.Time]; 

tstart = evtbl(strcmp(evtbl.EventName, 'image_1_adaptive'),:);
tend = evtbl(strcmp(evtbl.EventName, 'end_encoding'),:);
Tencode = [tstart.Time, tend.Time]; 

% art remove 
XArtRem = artremove(tRel, X, g, 15);

% filter, phase 
BPF = buildFIRBPF(fs,4,9);
Xf = filtfilt(BPF,1,XArtRem);
Xph = angle(hilbert(Xf));
figure('Theme','light'); 
plot(t, X-mean(XArtRem), 'k','LineWidth',2); hold on; grid on; 
plot(t, XArtRem-mean(XArtRem), 'b','LineWidth',1.75);
plot(t, Xf, 'm','LineWidth',1.5)
yrng = ylim();
yrng = 0.75*yrng;
yrng = [yrng; yrng]; yrng = yrng(:);
for ti = 1:height(Tencode)
    patch([Tencode(ti,:), fliplr(Tencode(ti,:))], yrng, ...
        'r', 'FaceAlpha',0.3, 'EdgeColor','none');
end
for ti = 1:height(Tdecode)
    patch([Tdecode(ti,:), fliplr(Tdecode(ti,:))], yrng, ...
        'g', 'FaceAlpha',0.3, 'EdgeColor','none');
end

% analysis 
phencode = []; phdecode = [];
for ti = 1:height(Tencode)
    [~,~,~,ph] = getPhaseEP(X,Xph,g,t, Tencode(ti,:),[],fs);
    phencode = [phencode; ph];
end
for ti = 1:height(Tdecode)
    [~,~,~,ph] = getPhaseEP(X,Xph,g,t, Tdecode(ti,:),[],fs);
    phdecode = [phdecode; ph];
end
figure; sgtitle('Stim Phase');
subplot(1,2,1); 
polarhistogram(phencode, 12);
phmean = circ_mean(phencode); phconf = circ_std(phencode);
phmean = phmean*180/pi; phconf = phconf*180/pi;
title('Encode'); subtitle([num2str(phmean),'±',num2str(phconf),'°']);
subplot(1,2,2); 
polarhistogram(phdecode, 12);
phmean = circ_mean(phdecode); phconf = circ_std(phdecode);
phmean = phmean*180/pi; phconf = phconf*180/pi;
title('Decode'); subtitle([num2str(phmean),'±',num2str(phconf),'°']);

%% helper(s)

function x = artremoveAR(mdl, mdlord, x, gidx, artstart, artend)
for gi = gidx'
    artidx = [artstart, artend] + gi;
    artidx(1) = max(1, artidx(1));
    artidx(2) = min(length(x), artidx(2));
    artidx = artidx(1):artidx(2);
    if numel(artidx)
        pastidx = artidx(1)-mdlord-1;
        if pastidx > 0
            pastidx = pastidx + (1:mdlord) - 1;
            xpast = x(pastidx); xpast = xpast-mean(xpast);
            x(artidx) = myFastForecastAR(mdl, xpast, length(artidx));
        end
    end
end
end

function XArtRem = artremove(tRel, X, g, artdur)
G = movmean(g,artdur)>eps;
Gshift = floor(artdur/2) - 1;
G = [false(Gshift,1); G(1:(end-Gshift))];
tRelNoArt = tRel(~G);
XNoArt = X(~G,:);
XArtRem = interp1(tRelNoArt, XNoArt, tRel, 'linear', 'extrap');
end

function getplotPhaseEP(x,xph,trig,t,trng,EPidx,fs)
[tEP,EPm,EPv, ph,phmean,phconf] = ...
    getPhaseEP(x,xph,trig,t,trng,EPidx,fs);
plotPhaseEP(tEP,EPm,EPv, ph,phmean,phconf);
end

function plotPhaseEP(tEP,EPm,EPv, ph,phmean,phconf)
phmean = phmean*180/pi; phconf = phconf*180/pi;
subplot(2,1,1); polarhistogram(ph, 12);
title('Stim Phase');
subtitle([num2str(phmean),'±',num2str(phconf),'°']);
subplot(2,1,2); 
errorbar(tEP,EPm,EPv); grid on;
xlabel('sec'); ylabel('uV'); title('Evoked Potential');
end

function [tEP,EPm,EPv, ph,phmean,phconf] = ...
    getPhaseEP(x,xph,trig,t,trng,EPidx,fs)
xidx = (t > trng(1)) & (t < trng(2));
x = x(xidx); xph = xph(xidx); trig = find(trig(xidx));
if ~isempty(EPidx)
    tEP = EPidx(1):EPidx(2); tEP = tEP/fs;
    EPl = diff(EPidx)+1;
    EP = nan(EPl,length(trig));
    for ti = 1:length(trig)
        EP(:,ti) = x(trig(ti)+(EPidx(1):EPidx(2)));
    end
    bldur = 5-EPidx(1);
    EPos = EP(1:bldur,:); EPos = mean(EPos);
    EP_ = EP - EPos; 
    EPm = mean(EP_,2); EPv = std(EP_,[],2);
else
    tEP = []; EPm = []; EPv = [];
end
ph = xph(trig); phmean = circ_mean(ph); phconf = circ_std(ph);
end