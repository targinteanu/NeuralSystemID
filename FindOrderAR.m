function [bestOrd, bestAIC, mdl] = ...
    FindOrderAR(dtaBL1ch, predHzn, fbnd, doResample, maxOrd, numTests)

if nargin < 6
    numTests = [];
    if nargin < 5
        maxOrd = [];
        if nargin < 4
            doResample = false;
            if nargin < 3
                fbnd = [];
                if nargin < 2
                    predHzn = [];
                end
            end
        end
    end
end

if isempty(maxOrd)
    maxOrd = 30;
end
if isempty(numTests)
    numTests = 100;
end
if isempty(predHzn)
    predHzn = .03; % s
end

try 
    Fs = dtaBL1ch.Properties.SampleRate;
catch ME
    warning(ME.message);
    Fs = nan;
end
if isnan(Fs) || isempty(Fs) 
    Fs = 1/mean(seconds(diff(dtaBL1ch.Time)));
end

%% filter, downsample 

FsNew = Fs;
if ~isempty(fbnd)
    tOrig = dtaBL1ch.Time;
    samplerat = 10; % replace with something based on nyquist rate
    bpf = buildFIRBPF(Fs, fbnd(1), fbnd(2));
    dtaBL1ch = FilterTimetable(@(f,x) filtfilt(f,1,x), bpf, dtaBL1ch);
    dtaOrig = dtaBL1ch;
    if doResample
        dtaBL1ch = downsampleTimetable(dtaBL1ch, samplerat);
        try 
            FsNew = dtaBL1ch.Properties.SampleRate;
        catch ME
            warning(ME.message);
            FsNew = nan;
        end
        if isnan(FsNew) || isempty(FsNew)
            FsNew = 1/mean(seconds(diff(dtaBL1ch.Time)));
        end
    end
end

%% compute

disp('Testing models on different training sets.')
progtick = .05; prog = 0;

L = floor(height(dtaBL1ch)/numTests);
datarat = L/maxOrd;
disp(['Data size is ',num2str(datarat),' times variable size.'])

l = ceil(L/2);
ts = linspace(l, height(dtaBL1ch)-l, numTests);
ts = round(ts);

W = nan(numTests, maxOrd+1);
r = 1;
for t = ts
    trng = t + [-1,1]*(l-1);
    y = dtaBL1ch(trng(1):trng(2),:);
    mdl = ar(y, maxOrd, 'ls');
    W(r,:) = mdl.A;
    r = r+1;
    prog = prog + 1/numTests;
    if prog > progtick
        disp([' - ',num2str(100*r/numTests),'% complete.'])
        prog = 0;
    end
end

%% plot
Wavg = mean(W);
Wstd = std(W);
Wsem = Wstd/sqrt(numTests);
Ts = tinv(.998, numTests-1);
Wp = zeros(size(Wavg));
for c = 1:width(W)
    [~,Wp(c)] = ttest(W(:,c));
end
figure; b = bar(Wavg);
hold on; grid on; 
errorbar(Wavg,Wsem*Ts,'.');
%{
text(b.XData, b.YData, string(Wp), ...
    "HorizontalAlignment","center", "VerticalAlignment","middle");
%}
xlabel('tap'); ylabel('coefficient');
legend('Expected', '99.9%CI');

%% compute AICs

disp('Testing models of different orders.')
progtick = .05; prog = 0;

AICtype = 'nAIC';

AIC = nan(1,maxOrd);
for ord = 1:maxOrd
    mdl = ar(y, ord, 'ls');
    AIC(ord) = aic(mdl, AICtype);
    prog = prog + 1/maxOrd;
    if prog > progtick
        disp([' - ',num2str(100*ord/maxOrd),'% complete.'])
        prog = 0;
    end
end

figure; bar(AIC); grid on;
xlabel('Model Order'); ylabel(AICtype); 
title('Model Information Criterion');

[bestAIC,bestOrd] = max(AIC);
mdl = ar(dtaBL1ch, bestOrd, 'yw');

%% forecast/predict

k = ceil(predHzn * FsNew);
LL = floor(height(dtaBL1ch)/2);
dtaFore = myFastForecastAR(mdl, dtaBL1ch{1:LL,1}, LL);
dtaFore = timetable(dtaBL1ch.Time(LL+(1:LL)), dtaFore);
dtaPred = myPredict(mdl, dtaBL1ch, k, true);

figure; 
plot(dtaOrig.Time, dtaOrig{:,1}, ':k', 'LineWidth', 1);
hold on; grid on;
plot(dtaBL1ch.Time, dtaBL1ch{:,1}, 'k'); 
plot(dtaPred.Time, dtaPred{:,1}, 'b');
plot(dtaFore.Time, dtaFore{:,1}, 'r');
[rho, p] = corr(dtaBL1ch{:,1}, dtaPred{:,1})

% upsample if applicable 
if (~isempty(fbnd)) && doResample
    dtaBL1ch = retime(dtaBL1ch,'regular','nearest','SampleRate',Fs);
    dtaPred  = retime(dtaPred, 'regular','nearest','SampleRate',Fs);
    dtaFore  = retime(dtaFore, 'regular','nearest','SampleRate',Fs);
    dtaBL1ch = FilterTimetable(@(f,x) filtfilt(f,1,x), bpf, dtaBL1ch);
    dtaPred  = FilterTimetable(@(f,x) filtfilt(f,1,x), bpf, dtaPred);
    dtaFore = FilterTimetable(@(f,x) filtfilt(f,1,x), bpf, dtaFore);

    plot(dtaBL1ch.Time, dtaBL1ch{:,1}, '--k'); hold on; grid on;
    plot(dtaPred.Time, dtaPred{:,1}, '--b');
    plot(dtaFore.Time, dtaFore{:,1}, '--r');
    [rho, p] = corr(dtaBL1ch{:,1}, dtaPred{:,1})
end
[rho, p] = corr(dtaOrig{:,1}, dtaPred{2:end,1})

end