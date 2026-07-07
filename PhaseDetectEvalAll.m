%% display params

mdlnames = {'Constant', 'Adaptive'};
bndnames = {'Beta', 'Theta'};
nbin = 48;
bedge = linspace(0,2*pi,nbin);
bedge = bedge - mean(bedge(1:2)); % center 0

clr = {[0.0660    0.4430    0.7450], ... blue 
       [0.8660    0.3290         0], ... red 
       [0.2310    0.6660    0.1960], ... green
       [0.5210    0.0860    0.8190], ... purple
       [0.6193    0.4627    0.0833], ... gold
       ...[0.9290    0.6940    0.1250], ... yellow
       [0.2588    0.5294    0.4471], ... teal
       [0.8190    0.0150    0.5450]  ... dark red
       };

%% load data from all files

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
    if ~((filedata.ARord==50) && (filedata.ARwin==1000) && ...
            (filedata.learnrate==0.05) && (filedata.packetSize)==20)
        warning('AR settings not consistent between files!')
    end
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

% ANOVA(-like) test : assign groups 
ERR3 = ERR(:,:,:,1); % ignore SD within subj/cond
GRP3 = (1:size(ERR3,1))' * ones(1,size(ERR3,2));
GRP3 = repmat(GRP3,[1,1,size(ERR3,3)]);

% calc errorbar values 
ERR2mv = nan([size(ERR3,1),size(ERR3,3),3]); 
for r = 1:size(ERR3,1)
    for c = 1:size(ERR3,3)
        ERR3rc = squeeze(ERR3(r,:,c))';
        N = length(ERR3(r,:,c));
        %ERR2mv(r,c,1) = mean(ERR3rc);
        %ERR2mv(r,c,2) = tinv(.995,N-1)*std(ERR3rc)/sqrt(N);
        ERR2mv(r,c,1) = circ_mean(ERR3rc);
        ERR2mv(r,c,2) = circ_confmean(ERR3rc, .05);
        ERR2mv(r,c,3) = circ_std(ERR3rc);
    end
end

%{
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
%}


% polar histogram of distributions with conf bars
figure('Position',[1 1 1125 825], 'WindowStyle','normal', ...
    'Theme','light', 'Color','w');
tiledlayout(1,size(ERR,3),'TileSpacing','compact');

for c = 1:size(ERR,3)
    ax2(c) = nexttile;
    ERR3c = ERR3(:,:,c); GRP3c = GRP3(:,:,c);
    R = zeros(1,size(ERR,1));
    for r = 1:size(ERR,1)
        phTgt = phTargets(r);
        ph = polarhistogram(ERR3c(r,:)+phTgt, 'BinEdges',bedge, 'FaceColor',clr{r}, ...
            'EdgeColor','none', 'FaceAlpha',0.4);
        hold on; %R = max(R, max(ph.Values));
        R(r) = max(ph.Values);
    end
    R = 1.1*R;
    for r = 1:size(ERR,1)
        phTgt = phTargets(r);
        polarboxplot(phTgt + ERR2mv(r,c,1), ERR2mv(r,c,2), ERR2mv(r,c,3), ...
            R(r), max(R), clr{r});
    end
    [p,ptbl] = circ_wwtest(ERR3c(:), GRP3c(:));

    ax2(c)=gca(); ax2(c).FontSize = 12;
    ax2(c).ThetaTick = 0:45:360;
    for ax2cti = 1:2:length(ax2(c).ThetaTickLabel)
        ax2(c).ThetaTickLabel{ax2cti} = '';
    end
    ax2(c).RTick = round(ax2(c).RTick(end)*[0.5,1]);
    title([mdlnames{c},' Model'], 'FontSize',16);
    subtitle(['Watson-Williams p value: ',num2str(p)], 'FontSize',14);
end
    lgd = "Target: "+string(phTargets*180/pi)+"°";
    %lgd = [lgd; repmat("mean ± 95% C.I.", size(lgd))];
    %legend(lgd(:), 'Location','northoutside');
    lgd = legend(lgd(:), 'FontSize',18);
    lgd.Layout.Tile = 'east';
sgtitle('Mean phase error by target phase', 'FontSize',20);


% histogram of squared error 
figure;
for c = 1:size(ERR,3)
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
figure('Position',[272 297 715 400], 'WindowStyle','normal', ...
    'Theme','light', 'Color','w');
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
selBnd = {selBeta, selTheta};
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

%% analysis by cond 

figure('Position',[1 1 1125 825], 'WindowStyle','normal', ...
    'Theme','light', 'Color','w');
tiledlayout(length(ERRbnd),size(ERR,3),'TileSpacing','compact');

for b = 1:length(ERRbnd)
    for m = 1:size(ERR,3)
        ax3(b,m) = nexttile;
        ERRbm = ERRbnd{b}(:,:,m);
        INFOb = INFObnd{b}; INFObCond = lower([INFOb.Cond]);

        selTask = (...
            contains(INFObCond, "task") | ...
            contains(INFObCond, "testing") | ...
            contains(INFObCond, "sternberg") | ...
            contains(INFObCond, "exam") | ...
            contains(INFObCond, "serialdigitalio") ) & ...
            ~contains(INFObCond, "baseline") ;
        selBaseline = strcmp(INFObCond, "baseline");
        selOther = (~selBaseline) & (~selTask);
        ERRbln = ERRbm(:,selBaseline); ERRbln = ERRbln(:);
        ERRtsk = ERRbm(:,selTask);     ERRtsk = ERRtsk(:);
        ERRoth = ERRbm(:,selOther);    ERRoth = ERRoth(:);

        p1a = circ_kuipertest(ERRbln, ERRtsk);
        p2a = circ_kuipertest(ERRbln, ERRoth);
        [~,p1b] = ttest2(ERRbln.^2, ERRtsk.^2, 'tail','left', 'Vartype','unequal');
        [~,p2b] = ttest2(ERRbln.^2, ERRoth.^2, 'tail','left', 'Vartype','unequal');

        %%{
        ERRbmmv = [circ_mean(ERRbln), circ_confmean(ERRbln), circ_std(ERRbln), rms(ERRbln); ...
                   circ_mean(ERRtsk), circ_confmean(ERRtsk), circ_std(ERRtsk), rms(ERRtsk); ...
                   circ_mean(ERRoth), circ_confmean(ERRoth), circ_std(ERRoth), rms(ERRoth)];

        ph = polarhistogram(ERRbln, 'BinEdges',bedge, 'FaceColor',clr{3}, ...
            'EdgeColor','none', 'FaceAlpha',0.4); 
        hold on; R = max(ph.Values);
        ph = polarhistogram(ERRtsk, 'BinEdges',bedge, 'FaceColor',clr{4}, ...
            'EdgeColor','none', 'FaceAlpha',0.4);
        R = max(R, max(ph.Values));
        ph = polarhistogram(ERRoth, 'BinEdges',bedge, 'FaceColor',clr{5}, ...
            'EdgeColor','none', 'FaceAlpha',0.4);
        R = max(R, max(ph.Values)); R = 1.25*R;
        %{
        polarregion(ERRbmmv(1,1) + [-1,1]*ERRbmmv(1,2), .5*[-1,1]+R, ...
            "FaceColor",clr{3}, "FaceAlpha",0.8, "EdgeColor",'k');
        polarregion(ERRbmmv(2,1) + [-1,1]*ERRbmmv(2,2), .5*[-1,1]+R, ...
            "FaceColor",clr{4}, "FaceAlpha",0.8, "EdgeColor",'k');
        polarregion(ERRbmmv(3,1) + [-1,1]*ERRbmmv(3,2), .5*[-1,1]+R, ...
            "FaceColor",clr{5}, "FaceAlpha",0.8, "EdgeColor",'k');
        %}
        for r = 1:3
            polarboxplot(ERRbmmv(r,1), ERRbmmv(r,2), ERRbmmv(r,3), ...
                R, R, clr{r+2});
        end

        ERRbmmv = ERRbmmv*180/pi;
        disp(['Model ',mdlnames{m},', Band ',bndnames{b},...
            ' - circmean, 95CI, std, rmse (deg): BL; task; other'])
        disp(ERRbmmv)
        %{
        lgd = [["Baseline"; "Task"; "Other"]; ...
            string(ERRbmmv(:,1))+"±"+string(ERRbmmv(:,2))]';
        legend(lgd(:), 'Location','eastoutside');
        %}
        %{
        legend(["Baseline"; "Task"; "Other"]...
            +", "+string(ERRbmmv(:,1))+"±"+string(ERRbmmv(:,2)), ...
            'Location','eastoutside');
        %}
        %}
        %{
        ERRbmmv = nan(3,2);
        [ERRbmmv(1,1), ERRbmmv(1,2)] = advPolHist(ERRbln, clr{3}, bedge);
        [ERRbmmv(2,1), ERRbmmv(2,2)] = advPolHist(ERRtsk, clr{4}, bedge);
        [ERRbmmv(3,1), ERRbmmv(3,2)] = advPolHist(ERRoth, clr{5}, bedge);
        ERRbmmv = ERRbmmv*180/pi;
        lgd = ["Baseline", "Task", "Other"]+" "+[""; "mean"; "95% C.I."];
        legend(lgd(:), 'Location','eastoutside');
        %}
        %{
        for c = 1:height(ERRbmmv)
            polarplot(ERRbmmv(c,1), 1, 'o', ...
                'LineWidth',4, 'Color',clr{c+2}, 'MarkerSize',10);
            hold on;
            polarplot(ERRbmmv(c,1)+[-1,1]*ERRbmmv(c,2), [1,1], ...
                'LineWidth',4, 'Color',clr{c+2});
        end
        lgd = ["Baseline", "Task", "Other"]+" "+["mean"; "95% C.I."];
        legend(lgd(:), 'Location','eastoutside');
        %}
        ax3(b,m) = gca(); ax3(b,m).FontSize = 12;
        ax3(b,m).ThetaTick = 0:45:360;
        for ax1bti = 1:2:length(ax3(b,m).ThetaTickLabel)
            ax3(b,m).ThetaTickLabel{ax1bti} = '';
        end
        ax3(b,m).RTick = round(ax3(b,m).RTick(end)*[0.5,1]);
        title([mdlnames{m},' Model, ',bndnames{b},' band'], 'FontSize',16);
        subtitle({['p1 = ',num2str(p1a),' | ',num2str(p1b)]; ...
                  ['p2 = ',num2str(p2a),' | ',num2str(p2b)]}, ...
                  'FontSize',14);
    end
end

sgtitle('Phase error by condition', 'FontSize',20);
lgd = legend({'Baseline', 'Task', 'Other'}, 'FontSize',18);
%lgd.Layout.Tile = 'eastoutside';

%% analysis of all 

figure('Position',[1 1 1125 425], 'WindowStyle','normal', ...
    'Theme','light', 'Color','w'); 
tiledlayout(1,length(bndnames)+1,'TileSpacing','compact');

% beta/theta/all
for b = 1:(length(bndnames)+1)
    ax1(b) = nexttile;
    if b > length(bndnames)
        bndname = 'Both'; ERRb = ERR;
    else
        bndname = bndnames{b}; ERRb = ERR(:,selBnd{b},:,:);
    end

    ERRbm = cell(length(mdlnames),1); R = 0;
    for m = 1:length(mdlnames)
        ERRbm_ = ERRb(:,:,m,1);
        ERRbm{m} = ERRbm_(:);
        ph = polarhistogram(ERRbm_, 'BinEdges',bedge, 'FaceColor',clr{m}, ...
            'EdgeColor','none', 'FaceAlpha',0.4);
        hold on;
        R = max(R, max(ph.Values));
    end
    %R = 1*R + 0.5*[-1,1];
    p1 = circ_kuipertest(ERRbm{1}, ERRbm{2});
    [~,p2] = ttest2(ERRbm{1}.^2, ERRbm{2}.^2, 'tail', 'right');

    ERRbstats = [cellfun(@circ_mean, ERRbm), ...
                 cellfun(@circ_confmean, ERRbm), ...
                 cellfun(@circ_std, ERRbm), ...
                 cellfun(@rms, ERRbm)];
    disp(['Band ',bndname,' circmean, 95CI, std, rmse (deg):'])
    disp(ERRbstats*180/pi);
    for m = 1:length(mdlnames)
        %{
        polarregion(ERRbstats(m,1) + [-1,1]*ERRbstats(m,2), R, ...
            "FaceColor",clr{m}, "FaceAlpha",0.8, "EdgeColor",'k');
        %}
        polarboxplot(ERRbstats(m,1), ERRbstats(m,2), ERRbstats(m,3), ...
            R, R, clr{m});
    end

    ax1(b)=gca(); ax1(b).FontSize = 12;
    ax1(b).ThetaTick = 0:45:360;
    for ax1bti = 1:2:length(ax1(b).ThetaTickLabel)
        ax1(b).ThetaTickLabel{ax1bti} = '';
    end
    ax1(b).RTick = round(ax1(b).RTick(end)*[0.5,1]);
    title([bndname,' band'], 'FontSize',16);
    subtitle(['p = ',num2str(p1),' | ',num2str(p2)], 'FontSize',14)
end

sgtitle('Phase Error', 'FontSize',20)
lgd = legend(string(mdlnames), 'FontSize',18);
lgd.Layout.Tile = 'east';

%% pie by cycle of extra/missing 

figure('Position',[1 1 875 625], 'WindowStyle','normal', ...
    'Theme','light', 'Color','w'); 
tiledlayout(length(mdlnames),length(bndnames)+1,'TileSpacing','compact');

for m = 1:length(mdlnames)
    NUMm = NUM(:,:,m,:);
    NUMm = sum(NUMm,1);
    
    % beta/theta
    for b = 1:length(bndnames)
        ax(m,b) = nexttile;
        NUMmb = NUMm(:,selBnd{b},:);
        NUMmb = sum(NUMmb,2);
        NUMmb = squeeze(NUMmb);
        hp = pie(ax(m,b), NUMmb);
        hpType = arrayfun(@(G) string(G.Type), hp);
        htxt = hp(strcmp(hpType, 'text'));
        for ht = htxt
            ht.FontSize = 14;
        end
        title([mdlnames{m},' model'], 'FontSize',16); 
        subtitle([bndnames{b},' band'], 'FontSize',16);
    end

    % all bnd
    ax(m,b+1) = nexttile;
    NUMm = sum(NUMm,2);
    NUMm = squeeze(NUMm);
    hp = pie(ax(m,b+1), NUMm);
    hpType = arrayfun(@(G) string(G.Type), hp);
    htxt = hp(strcmp(hpType, 'text'));
    for ht = htxt
        ht.FontSize = 14;
    end
    title([mdlnames{m},' model'], 'FontSize',16);
    subtitle('Both bands', 'FontSize',16);

end
    lgd = legend({'Missing', 'Extra', 'Correct'}, 'FontSize',18);
    lgd.Layout.Tile = 'east';

sgtitle('Number of Stimulations', 'FontSize',20)

%% helper(s) 

function polarboxplot(thetaC, thetaW1, thetaW2, rC, rMax, colr)
rW = rMax/20; 
% central polar box 
polarregion(thetaC + [-1,1]*thetaW1, rC + [-1,1]*rW, ...
    "FaceColor",colr, "EdgeColor",colr, "FaceAlpha",0.2, "LineWidth",2);
hold on;
% polar whiskers 
polarregion(thetaC + [-1,1]*thetaW2, [rC,rC], ...
    "FaceColor","None", "EdgeColor",colr, "LineWidth",2);
% whisker vertical indicators 
polarregion(thetaC - [1,1]*thetaW2, rC + [-1,1]*rW, ...
    "FaceColor","None", "EdgeColor",colr, "LineWidth",2);
polarregion(thetaC + [1,1]*thetaW2, rC + [-1,1]*rW, ...
    "FaceColor","None", "EdgeColor",colr, "LineWidth",2);
polarregion(thetaC + [0,0], rC + [-1,1]*rW, ...
    "FaceColor","None", "EdgeColor",colr, "LineWidth",2);
end

function [thetaAvg, thetaErb] = advPolHist(thetas, colr, binarg)

% Make a polar histogram of thetas that also indicates the circular mean
% and 95% CI above the histogram. 
% Optional second input dictates bin edges. If scalar, it will be
% considered #bins; otherwise, bin edges. 
% Optional third input determines color. 
if nargin < 3
    binarg = [];
    if nargin < 2
        colr = [];
    end
end
isbinedge = length(binarg) > 1;

% histogram 
if isempty(colr)
    if isempty(binarg)
        ph = polarhistogram(thetas);
    else
        if isbinedge
            ph = polarhistogram(thetas, "BinEdges",binarg);
        else
            ph = polarhistogram(thetas, binarg);
        end
    end
else
    if isempty(binarg)
        ph = polarhistogram(thetas, 'FaceColor',colr);
    else
        if isbinedge
            ph = polarhistogram(thetas, "BinEdges",binarg, 'FaceColor',colr);
        else
            ph = polarhistogram(thetas, binarg, 'FaceColor',colr);
        end
    end
end

% summary stats 
R = 1.25 * max(ph.Values);
thetaAvg = circ_mean(thetas);
%thetaErb = circ_confmean(thetas, .05);
thetaErb = circ_std(thetas);
if isempty(colr)
    colr = ph.FaceColor;
end
hold on; 
polarregion(thetaAvg + [-1,1]*thetaErb, [0,R], ...
    'FaceColor',colr)
%{
polarplot(thetaAvg, R, 'o', 'Color',colr, 'LineWidth',3, 'MarkerSize',8);
polarplot(thetaAvg + [-1,1]*thetaErb, [R,R], ...
    '-', 'Color',colr, 'LineWidth',3, 'MarkerSize',8);
%}

end



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