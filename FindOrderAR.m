function [bestOrd, bestAIC, bestMdl, rho, p] = ...
    FindOrderAR(dta, predHzn, fbnd, doResample, maxOrd, numTests)

% 
% Experiment with different AR model orders. 
% First, AR models of maximum order will be fit to different training sets
% selected from the input data, and the consistency between AR coefficients
% will be examined, with the expectation coefficients consistently close to
% zero are beyond the useful order of the model. 
% Next, AR models of increasing order will be fit to the same training set,
% and the AIC of each will be computed, with the expectation that the best
% order model will have lowest AIC. 
% Finally, signal forecast by the selected model will be compared with
% actual data as a measure of model accuracy. 
% 
% Inputs: 
%   dta: data to model from one channel without artifact as a timetable
%   predHzn: how much time (sec) ahead to evaluate model prediction
%            accuracy; default 0.03
%   fbnd: bandpass filter [low, high] frequency bounds; if empty or omitted, 
%         no filtering will be performed 
%   doResample: if true, data will be downsampled by 10 before fitting model
%   maxOrd: maximum AR order to test (default 30)
%   numTests: will attempt fitting on this many different training sets
%             selected from the data 
%   
% Outputs: 
%   bestOrd: best order defined by lowest AIC 
%   bestAIC: BIC value of order bestOrd model
%   bestMdl: an AR model of order bestOrd trained on all the data
%   rho: Pearson's correlation coefficient between actual data and data
%        predicted by the model predHzn seconds into the future. 
%   p: p value of the alternate hypothesis that rho is nonzero. 
% 

%% manage inputs 

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
    Fs = dta.Properties.SampleRate;
catch ME
    warning(ME.message);
    Fs = nan;
end
if isnan(Fs) || isempty(Fs) 
    Fs = 1/mean(seconds(diff(dta.Time)));
end

%% filter, downsample 

FsNew = Fs;
if ~isempty(fbnd)
    tOrig = dta.Time;
    samplerat = 10; % replace with something based on nyquist rate
    bpf = buildFIRBPF(Fs, fbnd(1), fbnd(2));
    dta = FilterTimetable(@(f,x) filtfilt(f,1,x), bpf, dta);
    dtaOrig = dta;
    if doResample
        dta = downsampleTimetable(dta, samplerat);
        try 
            FsNew = dta.Properties.SampleRate;
        catch ME
            warning(ME.message);
            FsNew = nan;
        end
        if isnan(FsNew) || isempty(FsNew)
            FsNew = 1/mean(seconds(diff(dta.Time)));
        end
    end
end

%% compute

disp('Testing models on different training sets.')
progtick = .05; prog = 0;

L = floor(height(dta)/numTests);
datarat = L/maxOrd;
disp(['Data size is ',num2str(datarat),' times variable size.'])

l = ceil(L/2);
ts = linspace(l, height(dta)-l, numTests);
ts = round(ts);

W = nan(numTests, maxOrd+1);
r = 1;
for t = ts
    trng = t + [-1,1]*(l-1);
    y = dta(trng(1):trng(2),:);
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

AICtype = 'BIC';

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

[bestAIC,bestOrd] = min(AIC);
mdl = ar(dta, bestOrd, 'yw');
bestMdl = mdl;

%% forecast/predict

k = ceil(predHzn * FsNew);
LL = floor(height(dta)/2);
dtaFore = myFastForecastAR(mdl, dta{1:LL,1}, LL);
dtaFore = timetable(dta.Time(LL+(1:LL)), dtaFore);
dtaPred = myPredict(mdl, dta, k, true);

figure; 
plot(dtaOrig.Time, dtaOrig{:,1}, ':k', 'LineWidth', 1);
hold on; grid on;
plot(dta.Time, dta{:,1}, 'k'); 
plot(dtaPred.Time, dtaPred{:,1}, 'b');
plot(dtaFore.Time, dtaFore{:,1}, 'r');
[rho, p] = corr(dta{:,1}, dtaPred{:,1});

% upsample if applicable 
if (~isempty(fbnd)) && doResample
    dta = retime(dta,'regular','nearest','SampleRate',Fs);
    dtaPred  = retime(dtaPred, 'regular','nearest','SampleRate',Fs);
    dtaFore  = retime(dtaFore, 'regular','nearest','SampleRate',Fs);
    dta = FilterTimetable(@(f,x) filtfilt(f,1,x), bpf, dta);
    dtaPred  = FilterTimetable(@(f,x) filtfilt(f,1,x), bpf, dtaPred);
    dtaFore = FilterTimetable(@(f,x) filtfilt(f,1,x), bpf, dtaFore);

    plot(dta.Time, dta{:,1}, '--k'); hold on; grid on;
    plot(dtaPred.Time, dtaPred{:,1}, '--b');
    plot(dtaFore.Time, dtaFore{:,1}, '--r');
    [rho, p] = corr(dta{:,1}, dtaPred{:,1});
end
[rho, p] = corr(dtaOrig{:,1}, dtaPred{2:end,1});

end