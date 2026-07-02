%% load data from all files
mdlnames = {'Constant', 'Dynamic'};

files = dir(fullfile('AdaptAR','*_PhaseDetect*.mat'));
DURDYN = []; DURCON = []; 
INFO = [];
ERR = []; 
% dim 1: target phase
% dim 2: subject/condition 
% dim 3: const vs dynamic AR model 
% dim 4: mean, SD
NUM = [];
% dim 1: target phase
% dim 2: subject/condition 
% dim 3: const vs dynamic AR model 
% dim 4: extra/missing/correct

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
        fErrCond_ = nan(size(fErr,1),1,size(fErr,3),2); 
        fNumCond_ = nan(size(fNum,1),1,size(fNum,3),length(fNum{1}));
        for m = 1:size(fErr,3)
            fErrCond_(:,:,m,1) = arrayfun(@(r)...
                circ_mean(cell2mat(fErrCond(r,:,m))'), ...
                1:height(fErr), 'UniformOutput',true);
            fErrCond_(:,:,m,2) = arrayfun(@(r)...
                circ_std(cell2mat(fErrCond(r,:,m))'), ...
                1:height(fErr), 'UniformOutput',true);
            fNumCond__ = arrayfun(@(r)...
                sum(cell2mat(fNum(r,:,m)')), ...
                1:height(fNum), 'UniformOutput',false);
            fNumCond_(:,:,m,:) = cell2mat(fNumCond__');
        end
        ERR = cat(2,ERR,fErrCond_);
        NUM = cat(2,NUM,fNumCond_);
    end

end

phTargets = filedata.phTargets;
%phTargets = phTargets*180/pi;

clear fi f filedata ifo ifoCond dims cond m
clear fErr fNum fErrCond fNumCond fErrCond_ fNumCond_ fNumCond__

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

clr = {[0.0660    0.4430    0.7450], ... blue 
       [0.8660    0.3290         0], ... red 
       [0.2310    0.6660    0.1960], ... green
       [0.5210    0.0860    0.8190], ... purple
       [0.6193    0.4627    0.0833], ... gold
       ...[0.9290    0.6940    0.1250], ... yellow
       [0.2588    0.5294    0.4471], ... teal
       [0.8190    0.0150    0.5450]  ... dark red
       };

% ANOVA(-like) test : assign groups 
ERR3 = ERR(:,:,:,1); % ignore SD within subj/cond
GRP3 = (1:size(ERR3,1))' * ones(1,size(ERR3,2));
GRP3 = repmat(GRP3,[1,1,size(ERR3,3)]);

% calc errorbar values 
ERR2mv = nan([size(ERR3,1),size(ERR3,3),2]); 
for r = 1:size(ERR3,1)
    for c = 1:size(ERR3,3)
        ERR3rc = squeeze(ERR3(r,:,c))';
        N = length(ERR3(r,:,c));
        %ERR2mv(r,c,1) = mean(ERR3rc);
        %ERR2mv(r,c,2) = tinv(.995,N-1)*std(ERR3rc)/sqrt(N);
        ERR2mv(r,c,1) = circ_mean(ERR3rc);
        ERR2mv(r,c,2) = circ_confmean(ERR3rc, .05);
        %ERR2mv(r,c,2) = circ_std(ERR3rc);
    end
end

% polar box-like plot 
figure;
for c = 1:size(ERR,3)
    subplot(1,size(ERR,3),c);
    ERR3c = ERR3(:,:,c); GRP3c = GRP3(:,:,c);
    for r = 1:size(ERR,1)
        phTgt = phTargets(r);
        polarplot(phTgt,1,'o', 'LineWidth',4, 'Color',clr{r}, 'MarkerSize',10);
        hold on; 
        polarplot(phTgt+ERR2mv(r,c,1),1,'s', 'LineWidth',4, 'Color',clr{r}, 'MarkerSize',10);
        polarplot(phTgt+ERR2mv(r,c,1)+[-1,1]*ERR2mv(r,c,2), [1,1], ...
            '-', 'LineWidth',3, 'Color',clr{r});
    end
    [p,ptbl] = circ_wwtest(ERR3c(:), GRP3c(:));
    title([mdlnames{c},': stim vs target phase']);
    subtitle(['Watson-Williams p value: ',num2str(p)]);
    legend('Target Phase for Stim', 'Mean Phase of Stim', '95% C.I.', ...
        'Location','northoutside');
end

% polar histogram of distributions 
nbin = 48;
bedge = linspace(0,2*pi,nbin);
bedge = bedge - mean(bedge(1:2)); % center 0
figure;
for c = 1:size(ERR,3)
    subplot(1,size(ERR,3),c);
    ERR3c = ERR3(:,:,c); GRP3c = GRP3(:,:,c);
    for r = 1:size(ERR,1)
        polarhistogram(ERR3c(r,:), 'BinEdges',bedge, 'FaceColor',clr{r});
        hold on;
    end
    legend("Target: "+string(phTargets*180/pi)+"°", ...
        'Location','northoutside');
    [p,ptbl] = circ_wwtest(ERR3c(:), GRP3c(:));
    title([mdlnames{c},' Model: mean phase error by target phase']);
    subtitle(['Watson-Williams p value: ',num2str(p)]);
end

% histogram of squared error 
figure;
for c = 1:width(ERR2)
    subplot(1,size(ERR,3),c);
    ERR3c = ERR3(:,:,c); GRP3c = GRP3(:,:,c);
    for r = 1:size(ERR,1)
        histogram(ERR3c(r,:).^2, 36, 'FaceColor',clr{r});
        hold on;
    end
    legend(string(phTargets*180/pi));
    p = wanova(ERR3c(:).^2, GRP3c(:));
    title([mdlnames{c},' Phase squared error by target phase']);
    subtitle(['ANOVA p value: ',num2str(p)]);
    xlabel('squared error (rad^2)'); ylabel('count');
end

ERR2mv = ERR2mv*180/pi;
% cartesian bar plot of mean+errorbar
figure;
b = bar(phTargets*180/pi,ERR2mv(:,:,1)', 'LineWidth',1);
set(gca, 'FontSize',12)
xlabel('Target phase (deg)', 'FontSize',14); 
ylabel('Stim Mean Phase Error (°)', 'FontSize',14)
title('Accuracy vs Target Phase', 'FontSize',18)
xticks(0:90:360)
hold on; grid on; 
for c = 1:length(b)
    errorbar(b(c).XEndPoints, b(c).YData, ...
        ERR2mv(:,c,2), ERR2mv(:,c,2), ...
        '.', 'Color',b(c).EdgeColor, 'LineWidth',2, 'CapSize',8);
end
legend('Constant Model', 'Dynamic Model', '95% C.I.', '95% C.I.',... '±1SD', '±1SD', ...
    'Location','northoutside', 'FontSize',14, 'Orientation','horizontal')

%% analysis by band 
selBeta = strcmp({INFO.Band}, "Beta");% & strcmp([INFO.Cond], "Baseline");
selTheta = strcmp({INFO.Band}, "Theta");% & strcmp([INFO.Cond], "Baseline");
ERRbnd  = {ERR(:,selBeta,:,1), ERR(:,selTheta,:,1)}; 
INFObnd = {INFO(selBeta),      INFO(selTheta)};
figure;
for m = 1:size(ERR,3)
    subplot(1,size(ERR,3),m);
    ERRb = cell(size(ERRbnd));
    for b = 1:length(ERRbnd)
        ERRbm = ERRbnd{b}(:,:,m);
        ERRb{b} = ERRbm(:);
        polarhistogram(ERRbm(:), 'BinEdges',bedge, 'FaceColor',clr{b+2}); 
        hold on;
    end
    p = circ_kuipertest(ERRb{1}, ERRb{2});
    title([mdlnames{m},' Model: mean phase error by band']);
    subtitle(['Kuiper p value: ',num2str(p)]);
    legend("Beta", "Theta", 'Location','northoutside');
end

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