%% display params

bndnames = {'Beta', 'Theta'};
mdlnames = {'Constant', 'Adaptive'};
nbin = 48;
bedge = linspace(0,2*pi,nbin);
bedge = bedge - mean(bedge(1:2)); % center 0

clr = {[0.0660    0.4430    0.7450], ... blue 
       [0.8660    0.3290         0], ... red 
       [0.2310    0.6660    0.1960], ... green
       [0.5210    0.0860    0.8190], ... purple
       ...[0.6193    0.4627    0.0833], ... gold
       [0.6110    0.4660    0.1250], ... gold
       ...[0.9290    0.6940    0.1250], ... yellow
       ...[0.2588    0.4471    0.5294], ... teal
       [     0    0.6390    0.6390], ... teal
       [0.8190    0.0150    0.5450], ... pink
       [0.3720    0.1050    0.0310], ... brown
       [0.7170    0.1920    0.1720], ... dark red
       [0.0070    0.3450    0.0540], ... dark green
       [0.0620    0.2580    0.5010], ... dark blue 
       [0.2422    0.1504    0.6603], ... dark indigo
       [0.0638    0.7446    0.7292], ... teal
       [0.9769    0.9839    0.0805]  ... yellow
       };
FaceAlpha = 0.6;

%% load data from all files

files = dir(fullfile('AdaptAR','NoiseAnalysis','*_PhaseDetect*.mat'));
INFO = [];
ERR = []; 
% dim 1: target phase
% dim 2: subject/condition 
% dim 3: noise level
% dim 4: const vs dynamic AR model 
% dim 5: mean, SD
NUM = [];
% dim 1: target phase
% dim 2: subject/condition 
% dim 3: noise level
% dim 4: const vs dynamic AR model 
% dim 5: extra/missing/correct

for fi = 1:length(files)
    f = files(fi);

    filedata = load(fullfile(f.folder,f.name));
    if ~((filedata.ARord==50) && (filedata.ARwin==1000) && ...
            (filedata.learnrate==0.05) && (filedata.packetSize)==20)
        warning('AR settings not consistent between files!')
    end

    ifo = [];
    ifo.Subj = f.name(1:8);
    if contains(f.name, 'Beta')
        ifo.Band = 'Beta';
    elseif contains(f.name, 'Theta')
        ifo.Band = 'Theta';
    else
        if upper(ifo.Subj(2)) == 'D'
            ifo.Band = 'Beta';
        elseif upper(ifo.Subj(2)) == 'Y'
            ifo.Band = 'Theta';
        else
            ifo.Band = '';
        end
    end
        ifo.Cond = "Baseline";
        fnOrig = split(f.name, '_');
        fnOrig = join(fnOrig(1:3), '_');
        fnOrig = [fnOrig{1},'.mat'];
        INFO = [INFO, ifo];

    fErr = filedata.errResults; fNum = filedata.numResults;
    dims = size(filedata.errResults);
    if length(dims) < 4
        dims(4)=1;
    end
    ERRm = []; NUMm = [];
    for m = 1:dims(4)
        % avg/sum across channels for this subj/cond
        fErrCond = fErr(:,:,:,m); fNumCond = fNum(:,:,:,m);
        fErrCond_ = nan(size(fErr,1),1,size(fErr,3),1,2); 
        fNumCond_ = nan(size(fNum,1),1,size(fNum,3),1,length(fNum{1}));
        for n = 1:size(fErr,3)
            fErrCond_(:,:,n,1,1) = arrayfun(@(r)...
                circ_mean(cell2mat(fErrCond(r,:,n))'), ...
                1:height(fErr), 'UniformOutput',true);
            fErrCond_(:,:,n,1,2) = arrayfun(@(r)...
                circ_std(cell2mat(fErrCond(r,:,n))'), ...
                1:height(fErr), 'UniformOutput',true);
            fNumCond__ = arrayfun(@(r)...
                sum(cell2mat(fNum(r,:,n,m)')), ...
                1:height(fNum), 'UniformOutput',false);
            fNumCond_(:,:,n,1,:) = cell2mat(fNumCond__');
        end
        ERRm = cat(4,ERRm,fErrCond_);
        NUMm = cat(4,NUMm,fNumCond_);
    end
    ERR = cat(2,ERR,ERRm); NUM = cat(2,NUM,NUMm);

end

noiselevels = filedata.noiselevels;
%noiselevels = 1./noiselevels.^2; % SNR 
phTargets = filedata.phTargets;
%phTargets = phTargets*180/pi;

clear fi f fnOrig filedata filedataOrig ifo ifoCond ifoDur dims m m
clear fErr fNum fErrCond fNumCond fErrCond_ fNumCond_ fNumCond__ ERRm NUMm

%% analysis by band 
selBeta = strcmp({INFO.Band}, "Beta");% & strcmp([INFO.Cond], "Baseline");
selTheta = strcmp({INFO.Band}, "Theta");% & strcmp([INFO.Cond], "Baseline");
selBnd = {selBeta, selTheta};
ERRbnd  = {ERR(:,selBeta,:,:), ERR(:,selTheta,:,:)}; 
INFObnd = {INFO(selBeta),      INFO(selTheta)};
%{
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
    %p = circ_kuipertest(ERRb{1}, ERRb{2});
    title({"Noise Amplitude (relative): "+string(noiselevels(m)); "mean phase error by band"});
    %subtitle(['Kuiper p value: ',num2str(p)]);
    legend("Beta", "Theta", 'Location','northoutside');
end
%}

%% aggregate/display all channel/target results 

errResultsAll_avg = nan(length(noiselevels),size(ERR,4));
errResultsAll_std = errResultsAll_avg;
for m = 1:size(ERR,4)
for n = 1:size(ERR,3)
    ERRn = ERR(:,:,n,m,1); ERRn = ERRn(:);
    errResultsAll_avg(n,m) = circ_mean(ERRn);
    errResultsAll_std(n,m) = circ_std(ERRn);
end
end
errResultsAll_avg = errResultsAll_avg*180/pi;
errResultsAll_std = errResultsAll_std*180/pi;

figure('Position',[272 297 715 400], 'WindowStyle','normal', ...
    'Theme','light', 'Color','w');
%bar(noiselevels, errResultsAll_avg); 
b = bar(errResultsAll_avg, 'LineWidth',2, 'EdgeColor','k', ...
    'FaceAlpha',FaceAlpha);%, 'FaceColor',clr{2}); 
for bi = 1:length(b)
    b(bi).FaceColor = clr{bi};
end
bX1 = {b.XData}; bX2 = {b.XEndPoints};
bX1 = cell2mat(bX1')'; bX2 = cell2mat(bX2')'; 
xticks(1:length(errResultsAll_avg)); xticklabels(string(noiselevels));
set(gca, 'FontSize',12)
hold on; grid on; 
%errorbar(noiselevels,errResultsAll_avg, errResultsAll_std,errResultsAll_std, '.');
errorbar(bX2, errResultsAll_avg, errResultsAll_std, ...
        '.', 'Color','k', 'LineWidth',2, 'CapSize',8);
xlabel('Noise Amplitude (relative)', 'FontSize',14); 
ylabel('Circ. Mean Phase Error (°)', 'FontSize',14); 
title('Accuracy vs Noise Level', 'FontSize',18)
legend([mdlnames, '±1 Circular S.D.'], ... 
    'Location','southoutside', 'FontSize',14, 'Orientation','horizontal');

%{

errResultsAll_avg = (ERR(:,:,:,1));
errResultsAll_std = (ERR(:,:,:,2));
errResultsAll_avg = mean(errResultsAll_avg,2);
errResultsAll_avg = mean(errResultsAll_avg,1);
errResultsAll_avg = squeeze(errResultsAll_avg)*180/pi;
errResultsAll_std = rms(errResultsAll_std,2);
errResultsAll_std = rms(errResultsAll_std,1);
errResultsAll_std = squeeze(errResultsAll_std)*180/pi;

figure('Position',[272 297 715 400], 'WindowStyle','normal', ...
    'Theme','light', 'Color','w');
%bar(noiselevels, errResultsAll_avg); 
bar(errResultsAll_avg, 'LineWidth',2, 'EdgeColor','k', ...
    'FaceAlpha',FaceAlpha, 'FaceColor',clr{2}); 
xticks(1:length(errResultsAll_avg)); xticklabels(string(noiselevels));
set(gca, 'FontSize',12)
hold on; grid on; 
%errorbar(noiselevels,errResultsAll_avg, errResultsAll_std,errResultsAll_std, '.');
errorbar(errResultsAll_avg, errResultsAll_std, ...
        '.', 'Color','k', 'LineWidth',2, 'CapSize',8);
xlabel('Noise Amplitude (relative)', 'FontSize',14); 
ylabel('Phase Error (°)', 'FontSize',14); 
title('Accuracy vs Learning Rate', 'FontSize',18)
legend('Circular Mean', '±1 Circular S.D.', ... 
    'Location','northoutside', 'FontSize',14, 'Orientation','horizontal');

%}

%% bar by cycle of extra/missing 

figure('Position',[1 1 875 625], 'WindowStyle','normal', ...
    'Theme','light', 'Color','w'); 
tiledlayout(length(bndnames)+1,length(mdlnames), "TileSpacing","compact")

% beta/theta
for b = 1:length(bndnames)
    NUMb = NUM(:,selBnd{b},:,:,:);
    NUMb = sum(NUMb,1); % across tgt phase
    NUMb = sum(NUMb,2); % across subject
    NUMb = squeeze(NUMb);
    NUMb = NUMb./sum(NUMb,3); % fraction of total cycles 
    for m = 1:size(NUMb,2)
    nexttile;
    NUMbm = squeeze(NUMb(:,m,:));
    hb = bar(100*NUMbm, 'LineWidth',2, 'FaceAlpha',FaceAlpha);
    for bi = 1:length(hb)
        hb(bi).FaceColor = clr{11+bi};
    end
    xticks(1:length(noiselevels)); xticklabels(string(noiselevels));
    ylabel('% of Cycles', 'FontSize',14); 
    title([mdlnames{m},' Model'], 'FontSize',16);
    subtitle([bndnames{b},' band'], 'FontSize',16);
    end
end

% all bnd
NUMb = NUM;
    NUMb = sum(NUMb,1);
    NUMb = sum(NUMb,2);
    NUMb = squeeze(NUMb);
    NUMb = NUMb./sum(NUMb,3);
    for m = 1:size(NUMb,2)
    nexttile;
    NUMbm = squeeze(NUMb(:,m,:));
    hb = bar(100*NUMbm, 'LineWidth',2, 'FaceAlpha',FaceAlpha);
    for bi = 1:length(hb)
        hb(bi).FaceColor = clr{11+bi};
    end
    xticks(1:length(noiselevels)); xticklabels(string(noiselevels));
    ylabel('% of Cycles', 'FontSize',14); 
xlabel('Noise Amplitude (relative)', 'FontSize',14); 
title([mdlnames{m},' Model'], 'FontSize',16);
subtitle('Both bands', 'FontSize',16);
    end

lgd = legend('Missing Stimulation', 'Stimulated Extra', 'Stimulated Correctly', ...
    'FontSize',14, 'Orientation','horizontal');
lgd.Layout.Tile = 'south';

sgtitle('Number of Stimulations', 'FontSize',20)

%% all bnd only

figure('Position',[272 297 875 400], 'WindowStyle','normal', ...
    'Theme','light', 'Color','w');

for m = 1:size(NUMb,2)
    nexttile;
    NUMbm = squeeze(NUMb(:,m,:));
    hb = bar(100*NUMbm, 'LineWidth',2, 'FaceAlpha',FaceAlpha);
    for bi = 1:length(hb)
        hb(bi).FaceColor = clr{11+bi};
    end
    xticks(1:length(noiselevels)); xticklabels(string(noiselevels));
    set(gca, 'FontSize',12)
    ylabel('% of Cycles', 'FontSize',14); 
xlabel('Noise Amplitude (relative)', 'FontSize',14); 
title([mdlnames{m},' Model'], 'FontSize',16);
grid on;
end

sgtitle('Number of Stimulations vs Noise Level', 'FontSize',18);

lgd = legend('Missing Stimulation', 'Stimulated Extra', 'Stimulated Correctly', ...
    'FontSize',14, 'Orientation','horizontal');%, 'location','southoutside');
lgd.Layout.Tile = 'south';

%% helper(s) 

function y = fixspacing(x)
minspacing = max(x)/5;
[z,xi] = sort(x, 'ascend'); % z = x(xi)
for zi = 2:length(z)
    if z(zi)-z(zi-1) < minspacing
        z(zi) = z(zi-1) + minspacing;
    end
end
zj = arrayfun(@(i) find(xi==i), 1:length(xi));
y = z(zj);
end

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