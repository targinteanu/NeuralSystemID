DURDYN = []; DURCON = [];
ERR = {}; NUM = {}; TRG = {};
% dim 1: target phase
% dim 2: channel 
% dim 3: const vs dynamic AR model 

%% load data from all files

files = dir(fullfile('AdaptAR','*_PhaseDetect*.mat'));
for f = files'
    filedata = load(fullfile(f.folder,f.name));
    DURDYN = [DURDYN,filedata.durDyn];
    DURCON = [DURCON,filedata.durConst];
    ERR = cat(2,ERR,filedata.errResults);
    NUM = cat(2,NUM,filedata.numResults);
    TRG = cat(2,TRG,filedata.trgResults);
end

phTargets = filedata.phTargets;
phTargets = phTargets*180/pi;

%% computation time 
figure; 
histogram(DURCON); hold on; grid on; histogram(DURDYN);
set(gca, 'FontSize',12)
xlabel('duration (s)', 'FontSize',14); ylabel('count', 'FontSize',14);
title('Computation Time (single loop)', 'FontSize',18);
legend('Constant', 'Dynamic', 'location','northoutside', ...
    'Orientation','horizontal', 'FontSize',14);

%% cat across chans/subjs

ERR2 = cell(size(ERR,1), size(ERR,3));
NUM2 = cell(size(NUM,1), size(NUM,3));

for r = 1:size(ERR,1)
    for c = 1:size(ERR,3)
        ERR2{r,c} = cell2mat(ERR(r,:,c));
        NUM2{r,c} = cell2mat(NUM(r,:,c)');
    end
end

ERR2mv = nan([size(ERR2),2]); 
%NUM2mv = nan([size(NUM2),2]);
for r = 1:size(ERR2,1)
    for c = 1:size(ERR2,2)
        N = length(ERR2{r,c});
        ERR2mv(r,c,1) = mean(ERR2{r,c}*180/pi);
        ERR2mv(r,c,2) = tinv(.995,N-1)*std(ERR2{r,c}*180/pi)/sqrt(N);
        %N = length(NUM2{r,c});
        %NUM2mv(r,c,1) = mean(NUM2{r,c});
        %NUM2mv(r,c,2) = tinv(.995,N-1)*std(NUM2{r,c})/sqrt(N);
    end
end

figure;
b = bar(phTargets,ERR2mv(:,:,1)', 'LineWidth',1);
set(gca, 'FontSize',12)
xlabel('Target phase (deg)', 'FontSize',14); 
ylabel('Stim Phase Error (causal - offline) (deg)', 'FontSize',14)
title('Accuracy vs Target Phase', 'FontSize',18)
xticks(0:90:360)
hold on; grid on; 
for c = 1:length(b)
    errorbar(b(c).XEndPoints, b(c).YData, ...
        ERR2mv(:,c,2), ERR2mv(:,c,2), ...
        '.', 'Color',b(c).EdgeColor, 'LineWidth',2, 'CapSize',8);
end
legend('Constant', 'Dynamic', '99% C.I.', '99% C.I.',... '±1SD', '±1SD', ...
    'Location','northoutside', 'FontSize',14, 'Orientation','horizontal')

%% cat across tgt phases 

ERR3 = cell(1, size(ERR,3));
NUM3 = cell(1, size(NUM,3));

for c = 1:size(ERR,3)
    ERR3{1,c} = cell2mat(ERR2(:,c)');
    NUM3{1,c} = sum(cell2mat(NUM2(:,c)));
end

figure; 
polarhistogram(ERR3{1}, 256, 'EdgeColor','none'); 
hold on; polarhistogram(ERR3{2}, 256, 'EdgeColor','none');
title('Stim error (causal - offline)'); 
legend( ...
    ['Constant - RMSE ',num2str(rms(ERR3{1}))], ... 
    ['Dynamic - RMSE ',num2str(rms(ERR3{2}))], ...
    'Location','northoutside');
[~,p] = ttest2(ERR3{1}.^2, ERR3{2}.^2, 'tail', 'right');
subtitle(['p = ',num2str(p)])

figure; 
tiledlayout(1,2,'TileSpacing','compact');
ax1 = nexttile;
pie(ax1, NUM3{1});
title('Num. Stim. - Constant');
ax2 = nexttile;
pie(ax2, NUM3{2});
title('Num. Stim. - Dynamic');
lgd = legend({'Missing', 'Extra', 'Correct'});
lgd.Layout.Tile = 'east';