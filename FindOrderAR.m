function FindOrderAR(dtaBL1ch, fbnd, maxOrd, numTests)

if nargin < 4
    numTests = [];
    if nargin < 3
        maxOrd = [];
        if nargin < 2
            fbnd = [];
        end
    end
end

if isempty(maxOrd)
    maxOrd = 30;
end
if isempty(numTests)
    numTests = 100;
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

if ~isempty(fbnd)

tOrig = dtaBL1ch.Time;
samplerat = 10; % replace with something based on nyquist rate
bpf = buildFIRBPF(Fs, fbnd(1), fbnd(2));
dtaBL1ch = FilterTimetable(@(f,x) filtfilt(f,1,x), bpf, dtaBL1ch);
dtaOrig = dtaBL1ch;
dtaBL1ch = downsampleTimetable(dtaBL1ch, samplerat);

end

%% compute
L = floor(height(dtaBL1ch)/numTests);
datarat = L/maxOrd;
disp(['Data size is ',num2str(datarat),' times variable size.'])

l = ceil(L/2);
ts = linspace(l, height(dtaBL1ch)-l, numTests);
ts = round(ts);

W = nan(numTests, maxOrd+1);
r = 1;
for t = ts
    %tic
    trng = t + [-1,1]*(l-1);
    y = dtaBL1ch(trng(1):trng(2),:);
    mdl = ar(y, maxOrd, 'ls');
    W(r,:) = mdl.A;
    r = r+1;
    %toc
    %pause(.001)
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

%% forecast/predict

LL = floor(height(dtaBL1ch)/2);
dtaFore = myFastForecastAR(mdl, dtaBL1ch{1:LL,1}, LL);
dtaFore = timetable(dtaBL1ch.Time(LL+(1:LL)), dtaFore);
dtaPred = myPredict(mdl, dtaBL1ch, 3, true);

figure; 
plot(dtaOrig.Time, dtaOrig{:,1}, ':k', 'LineWidth', 1);
hold on; grid on;
plot(dtaBL1ch.Time, dtaBL1ch{:,1}, 'k'); 
plot(dtaPred.Time, dtaPred{:,1}, 'b');
plot(dtaFore.Time, dtaFore{:,1}, 'r');
[rho, p] = corr(dtaBL1ch{:,1}, dtaPred{:,1})

% upsample if applicable 
if ~isempty(fbnd)
    dtaBL1ch = retime(dtaBL1ch,'regular','nearest','SampleRate',Fs);
    dtaPred  = retime(dtaPred, 'regular','nearest','SampleRate',Fs);
    dtaFore  = retime(dtaFore, 'regular','nearest','SampleRate',Fs);
    dtaBL1ch = FilterTimetable(@(f,x) filtfilt(f,1,x), bpf, dtaBL1ch);
    dtaPred  = FilterTimetable(@(f,x) filtfilt(f,1,x), bpf, dtaPred);
    dtaFore = FilterTimetable(@(f,x) filtfilt(f,1,x), bpf, dtaFore);
end

plot(dtaBL1ch.Time, dtaBL1ch{:,1}, '--k'); hold on; grid on;
plot(dtaPred.Time, dtaPred{:,1}, '--b');
plot(dtaFore.Time, dtaFore{:,1}, '--r');
[rho, p] = corr(dtaBL1ch{:,1}, dtaPred{:,1})
[rho, p] = corr(dtaOrig{:,1}, dtaPred{2:end,1})

end