%% load data from all files

files = dir(fullfile('AdaptAR','*_PhaseDetect*.mat'));
DURDYN = []; DURCON = []; 
INFO = [];
ERR = {}; NUM = {}; 
% dim 1: target phase
% dim 2: subject/condition 
% dim 3: const vs dynamic AR model 

for fi = 1:length(files)
    f = files(fi);

    filedata = load(fullfile(f.folder,f.name));
    DURDYN = [DURDYN,filedata.durDyn];
    DURCON = [DURCON,filedata.durConst];

    ifo = [];
    ifo.Subj = f.name(1:8);
    if contains(f.name, 'Beta')
        ifo.Band = 'Beta';
    elseif contains(f.name, 'Theta')
        ifo.Band = 'Theta';
    else
        ifo.Band = '';
    end
    if contains(f.name, 'Nonbaseline')
        ifoCond = filedata.tblsToTest2Descs;
        dims = length(ifoCond);
        ifoCond = split(ifoCond, ":");
        if dims == 1
            ifoCond = ifoCond(1);
        else
            ifoCond = ifoCond(:,1);
        end
    else
        ifoCond = "Baseline";
    end

    fErr = filedata.errResults; fNum = filedata.numResults;
    dims = size(filedata.errResults);
    if length(dims) < 4
        dims(4)=1;
    end
    for cond = 1:dims(4)
        ifo.Cond = ifoCond(cond);
        INFO = [INFO, ifo];
        % avg/sum across channels for this subj/cond
        fErrCond = fErr(:,:,:,cond); fNumCond = fNum(:,:,:,cond);
        fErrCond_ = cell(size(fErr,1),1,size(fErr,3)); fNumCond_ = fErrCond_;
        for m = 1:size(fErr,3)
            fErrCond_(:,:,m) = arrayfun(@(r)...
                circ_mean(cell2mat(fErrCond(r,:,m))'), ...
                1:height(fErr), 'UniformOutput',false);
            fNumCond_(:,:,m) = arrayfun(@(r)...
                sum(cell2mat(fNum(r,:,m)')), ...
                1:height(fNum), 'UniformOutput',false);
        end
        ERR = cat(2,ERR,fErrCond_);
        NUM = cat(2,NUM,fNumCond_);
    end

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

disp(['Median time constant: ',num2str(median(DURCON))])
disp(['Median time dynamic: ',num2str(median(DURDYN))])

%% analysis by tgt phase 

% cat across chans/subjs
ERR2 = cell(size(ERR,1), size(ERR,3));
NUM2 = cell(size(NUM,1), size(NUM,3));
for r = 1:size(ERR,1)
    for c = 1:size(ERR,3)
        ERR2{r,c} = cell2mat(ERR(r,:,c));
        NUM2{r,c} = cell2mat(NUM(r,:,c)');
    end
end

% ANOVA test 
% cat across tgt phases 
ERR3 = cell(1, size(ERR,3));
NUM3 = cell(1, size(NUM,3));
GRP3 = cell(1, size(ERR,3));
for c = 1:size(ERR,3)
    ERR3{1,c} = cell2mat(ERR2(:,c)');
    NUM3{1,c} = sum(cell2mat(NUM2(:,c)));
    for r = 1:height(ERR2)
        GRP3{1,c} = [GRP3{1,c}, r*ones(size(ERR2{r,c}))];
    end
end

% kuiper test 
P = nan([height(ERR2),height(ERR2),width(ERR2)]);
for c = 1:width(ERR2)
    for r = 1:height(ERR2)
        for rr = (r+1):height(ERR2)
            P(r,rr,c) = circ_kuipertest(ERR2{r,c}',ERR2{rr,c}');
            %[~,P(r,rr,c)] = circ_ztest(ERR2{r,c}',ERR2{rr,c}');
            %[~,P(r,rr,c)] = kstest2((ERR2{r,c}.^2)', (ERR2{rr,c}.^2)');
            %[~,P(r,rr,c)] = ttest2((ERR2{r,c}.^2)', (ERR2{rr,c}.^2)', 'VarType','unequal');
        end
    end
end
P

% polar histogram of distributions 
figure;
for c = 1:width(ERR2)
    subplot(1,width(ERR2),c);
    for r = 1:height(ERR2)
        polarhistogram(ERR2{r,c},36);
        hold on;
    end
    legend(string(phTargets));
    [p,ptbl] = circ_wwtest(ERR3{1,c}', GRP3{1,c}');
    title(['ww p value: ',num2str(p)]);
end

% histogram of squared error 
figure;
for c = 1:width(ERR2)
    subplot(1,width(ERR2),c);
    for r = 1:height(ERR2)
        histogram(ERR2{r,c}.^2,36);
        hold on;
    end
    legend(string(phTargets));
    p = wanova((ERR3{1,c}.^2)', GRP3{1,c}');
    title(['Anova p value: ',num2str(p)]);
end

% calc errorbar values (CIs)
ERR2mv = nan([size(ERR2),2]); 
%NUM2mv = nan([size(NUM2),2]);
for r = 1:size(ERR2,1)
    for c = 1:size(ERR2,2)
        N = length(ERR2{r,c});
        %ERR2mv(r,c,1) = mean(ERR2{r,c}*180/pi);
        %ERR2mv(r,c,2) = tinv(.995,N-1)*std(ERR2{r,c}*180/pi)/sqrt(N);
        ERR2mv(r,c,1) = circ_mean(ERR2{r,c}');
        ERR2mv(r,c,2) = circ_confmean(ERR2{r,c}', .05);
        %ERR2mv(r,c,2) = circ_std(ERR2{r,c}');
        %N = length(NUM2{r,c});
        %NUM2mv(r,c,1) = mean(NUM2{r,c});
        %NUM2mv(r,c,2) = tinv(.995,N-1)*std(NUM2{r,c})/sqrt(N);
    end
end
ERR2mv = ERR2mv*180/pi;

% cartesian bar plot of mean+errorbar
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
legend('Constant', 'Dynamic', '95% C.I.', '95% C.I.',... '±1SD', '±1SD', ...
    'Location','northoutside', 'FontSize',14, 'Orientation','horizontal')

%% analysis of all 

% polar histo with p value 
figure; 
polarhistogram(ERR3{1}, 36);% 256, 'EdgeColor','none'); 
hold on; polarhistogram(ERR3{2}, 36);% 256, 'EdgeColor','none');
title('Stim error (causal - offline)'); 
legend( ...
    ['Constant - RMSE ',num2str(rms(ERR3{1}))], ... 
    ['Dynamic - RMSE ',num2str(rms(ERR3{2}))], ...
    'Location','northoutside');
%[~,p] = ttest2(ERR3{1}.^2, ERR3{2}.^2, 'tail', 'right');
%[~,p] = circ_ztest(ERR3{1}', ERR3{2}');
p = circ_kuipertest(ERR3{1}', ERR3{2}');
subtitle(['p = ',num2str(p)])

% pie by cycle of extra/missing 
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

%%
[~,p] = ttest2(ERR3{1}.^2, ERR3{2}.^2, 'tail', 'right')
[~,p] = circ_ztest(ERR3{1}', ERR3{2}')
[circ_mean(ERR3{1}'), mean(ERR3{1}'), median(ERR3{1}'), circ_zconf(ERR3{1}',.05), circ_confmean(ERR3{1}',.05)] *180/pi
[circ_mean(ERR3{2}'), mean(ERR3{2}'), median(ERR3{2}'), circ_zconf(ERR3{2}',.05), circ_confmean(ERR3{2}',.05)] *180/pi

%% helper(s) 

function CI = circ_zconf(x, p)
Z = norminv(1-(p/2));
N = length(x);
sx = sin(x); cx = cos(x);
R = sqrt(sum(cx)^2 + sum(sx)^2);
Rbar = R/N; % could use circ_r(x) instead 
CI = Z / ( Rbar * sqrt(2*N) );
end

function [z_stat, p_value] = circ_ztest(sample1, sample2)
    % Inputs: sample1, sample2 - Vectors of angles in radians (can be unequal length)
    
    N1 = length(sample1);
    N2 = length(sample2);
    
    % --- Group 1 Descriptive Stats ---
    sum_cos1 = sum(cos(sample1));
    sum_sin1 = sum(sin(sample1));
    theta_bar1 = atan2(sum_sin1, sum_cos1);
    R_bar1 = sqrt(sum_cos1^2 + sum_sin1^2) / N1;
    
    % --- Group 2 Descriptive Stats ---
    sum_cos2 = sum(cos(sample2));
    sum_sin2 = sum(sin(sample2));
    theta_bar2 = atan2(sum_sin2, sum_cos2);
    R_bar2 = sqrt(sum_cos2^2 + sum_sin2^2) / N2;
    
    % --- Directional Standard Errors (Unpooled) ---
    se1_sq = 1 / (2 * N1 * (R_bar1^2));
    se2_sq = 1 / (2 * N2 * (R_bar2^2));
    
    % --- Calculate Z-Statistic ---
    angular_diff = theta_bar1 - theta_bar2;
    z_stat = sin(angular_diff) / sqrt(se1_sq + se2_sq);
    
    % --- Two-Tailed p-value ---
    p_value = 2 * (1 - normcdf(abs(z_stat)));
    
end