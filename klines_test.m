x1 = randn(10,1)*2 + 5;
x2 = randn(12,1)*3 - 1;
x3 = randn(8,1)*.5 + 2;

y1 = 5*x1 + 6;
y2 = -.5*x2 - 2;
y3 = 4*x3 + 2;

y1 = y1+2*rand(size(y1));
y2 = y2+2*rand(size(y2));
y3 = y3+2*rand(size(y3));

colr = {'r','g','b'};

figure;
plot(x1,y1, ['.',colr{1}]); grid on; hold on;
plot(x2,y2, ['.',colr{2}]);
plot(x3,y3, ['.',colr{3}]);

XY = [x1,y1; x2,y2; x3,y3];
[labels, lines] = klines(XY, 3, 10000);

gm = fitgmdist(XY,3,...
    'CovarianceType','full',...
    'Replicates',20);

labels2 = cluster(gm,XY);

for lbl = 1:3
    xy = XY(labels==lbl,:);
    plot(xy(:,1), xy(:,2), ['o',colr{lbl}]);
    xy = XY(labels2==lbl,:);
    plot(xy(:,1), xy(:,2), ['d',colr{lbl}]);
end