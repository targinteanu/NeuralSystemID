maxOrd = 20;
numTests = 10;
L = floor(height(dtaBL1ch)/numTests);
datarat = L/maxOrd;
disp(['Data size is ',num2str(datarat),' times variable size.'])

l = ceil(L/2);
ts = linspace(l, height(dtaBL1ch)-l, numTests);
ts = round(ts);

W = nan(numTests, maxOrd+1);
r = 1;
for t = ts
    tic
    trng = t + [-1,1]*(l-1);
    y = dtaBL1ch(trng(1):trng(2),:);
    mdl = ar(y, maxOrd, 'ls');
    W(r,:) = mdl.A;
    r = r+1;
    toc
    pause(.001)
end

Wavg = mean(W);
Wstd = std(W);
figure; b = bar(Wavg);
hold on; grid on; 
errorbar(Wavg,3*Wstd,'.');
xlabel('tap'); ylabel('coefficient');
legend('Expected', '3SD');