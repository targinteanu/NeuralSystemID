maxOrd = 50;
numTests = 100;
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
    pause(.001)
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

%% forecast 
LL = floor(height(dtaBL1ch)/2);
dtaPred = myFastForecastAR(mdl, dtaBL1ch{1:LL,1}, LL);
figure; plot(dtaBL1ch{:,1}); grid on; hold on; 
plot([dtaBL1ch{1:LL,1}; dtaPred]);