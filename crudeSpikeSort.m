%% threshold definition 

% identify threshold(s) 
% alternatively do this based on mean/SD
[OL, lTH, uTH, mid] = isoutlier(x); 
OL = find(OL);
lOL = OL( x(OL) < mid );
uOL = OL( x(OL) > mid );

% decide whether to use upper or lower threshold 
% based on count or magnitude 
sgn = false; % flip signal sign
%{
if length(uOL) > length(lOL)
    TH = uTH;
elseif length(lOL) > length(uOL)
    TH = lTH;
    sgn = true;
else
    lD = mean(x(lOL)) - mid;
    uD = mean(x(uOL)) - mid;
    if abs(uD) > abs(lD)
        TH = uTH;
    else
        TH = lTH;
        sgn = true;
    end
end
%}
%[~,~,~,Tstat] = ttest2(abs(x(uOL)), abs(x(lOL)), 'Vartype','unequal');
%if Tstat.tstat > 0
xu = x(uOL).^2; xl = x(lOL).^2;
if mean(xu)/std(xu) > mean(xl)/std(xl)
    TH = uTH; % Set threshold to upper threshold
else
    TH = lTH; % Set threshold to lower threshold
    sgn = true; % Flip signal sign
end

if sgn
    TH = -TH;
    x = -x;
end

%% spike detection 
% use findpeaks 
% alternatively, look for positive thresh crossings 
[spkVal, spkIdx] = findpeaks(x); 
% consider adding minPeakProminence/Distance related to WFlen
spkSel = spkVal > TH;
spkIdx = spkIdx(spkSel); % Select indices of detected spikes
spkVal = spkVal(spkSel); % Select values of detected spikes

if sgn
    TH = -TH;
    x = -x;
    spkVal = -spkVal;
end

%% build waveform list 
WFdur = 0.004; % full duration, seconds
fs = 1/median(diff(t)); % hz
WFlen = ceil(WFdur/2 * fs); % half-duration, samples 
WF = nan(length(spkIdx), 2*WFlen+1);
for n = 1:length(spkIdx)
    idx = spkIdx(n);
    if (idx > WFlen) && (idx+WFlen <= length(x))
        WF(n,:) = x((idx-WFlen):(idx+WFlen));
    end
end
inan = sum(isnan(WF),2);
WF = WF(~inan,:);
tWF = (-WFlen:WFlen)/fs;

%% reduce dimension 
maxPC = 200; % max # components to consider 
thrE = 95; % cumulative %var explained target 
[C,WFS,L,T,E,M] = pca(WF);
EE = cumsum(E);
nPC = find(EE > thrE);
if isempty(nPC)
    nPC = maxPC;
else
    nPC = min(nPC);
    nPC = min(nPC, maxPC);
end
figure; plot(EE); hold on; grid on; plot(nPC, EE(nPC), 'o');
xlabel('# PCs'); ylabel('Cumulative % var explained');
title('PCA dim reduction');
nPC = min(nPC, size(WF,2));
WFS = WFS(:,1:nPC);

%% cluster waveforms 

% identify # clusters, k
krange = 2:ceil(size(WFS,1)/100);
%{
% Determine the optimal number of clusters using the elbow method
distances = pdist(WFS, 'euclidean');
linkageTree = linkage(distances, 'ward');
clusterDendrogram = dendrogram(linkageTree, 0);
%}
% Calculate the optimal number of clusters using the silhouette method
%silhouetteValues = silhouette(WFS, clusterIdx);
%optimalK = evalclusters(WFS, 'gmdistribution', 'Silhouette', 'KList', krange);
%k = optimalK.OptimalK;
%figure; plot(optimalK.InspectedK, optimalK.CriterionValues);
CriterionValues = nan(size(krange));
parfor ki = 57:length(krange)
    %tic
    [kidx,kc,kd] = kgmm(WFS, krange(ki));
    CriterionValues(ki) = mean(silhouette(WFS, kidx));
    %toc
end
figure; plot(krange, CriterionValues);
[~,ki] = max(CriterionValues); k = krange(ki);
grid on; hold on; plot(k, CriterionValues(krange==k), 'o');
xlabel('# clusters (k)'); ylabel('mean Silhouette score');
%{
dd = [];
ddd = pdist(WFS, 'euclidean'); ddd = mean(ddd);
for k = krange
    [ki,kc,kd] = kmeans(WFS, k);
    kn = arrayfun(@(kk) sum(ki==kk), 1:k)';
    kd = kd./kn; kd = mean(kd);
    dd = [dd, kd];
end
dd = dd/dd(1);
figure; plot(krange, dd);
%}

%% visualize clusters 
%[kidx,kc,kd] = kmeans(WFS, k);
[kidx,kc,kd] = kgmm(WFS, k);

% scatter 
figure;
for ki = unique(kidx)'
    WFi = WFS(kidx == ki, 1:3);
    plot3(WFi(:,1), WFi(:,2), WFi(:,3), '.'); hold on;
end
grid on;
title('Cluster scatterplot in PC space'); 
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

% waveform 
figure; 
errorbar(tWF, mean(WF), std(WF), 'Color',[.5 .5 .5]); hold on;
for ki = unique(kidx)'
    WFi = WF(kidx == ki, :);
    errorbar(tWF, mean(WFi), std(WFi)); hold on;
end
grid on;
title('Cluster Waveforms');
xlabel('time (s)'); ylabel('signal +- 1SD');

% timing
figure; 
plot(t, x); hold on; 
for ki = unique(kidx)'
    spki = spkIdx(kidx == ki);
    plot(t(spki), x(spki), 'o');
end
plot([t(1), t(end)], [TH, TH], '--');
grid on; 
xlabel('time (s)'); 

% ISI histo
figure; 
histogram(diff(spkIdx)/fs, 'FaceColor',[.5 .5 .5]); hold on;
for ki = unique(kidx)'
    spki = spkIdx(kidx == ki);
    histogram(diff(spki)/fs);
end
grid on;
title('Inter-Spike Interval');
xlabel('duration (s)'); ylabel('count');

%% auto/cross corr
figure; 
p = 1;
for ki = unique(kidx)'
    spki = spkIdx(kidx == ki);
    for kj = unique(kidx)'
        spkj = spkIdx(kidx == kj);
        subplot(k,k,p);
        [bincent, binvals] = corgm(spki, spkj, ceil(0.2*fs));
        bincent = bincent/fs; binvals = binvals*fs;
        bar(bincent, binvals);
        xlabel('delay (s)'); ylabel('avg # per s');
        title([num2str(kj),' given ',num2str(ki)]);
        p = p+1;
    end
end

%% helpers 


function [idx, C, S] = kgmm(X, k)
% fit Gaussian mixture model to WFS to obtain cluster assignments kidx
%rng(0); % for reproducibility
options = statset('MaxIter',10000);

% try to fit GMM with k components, using regularization to avoid singular covariances
gm = fitgmdist(X, k, 'Options', options, 'RegularizationValue', 1e-6, ...
    'Replicates', 5, 'CovarianceType', 'full');

% cluster indices: assign each point to the component with highest posterior
idx = cluster(gm, X);

% also compute cluster centers in PC space (optional for plotting/inspection)
C = gm.mu;
S = gm.Sigma;
end


function [bincent, binvals] = corgm(y, z, trng)
% compute cross-event conditional occurrence (histogram of z around y)
% y and z are vectors of timepoints (in same units as trng), trng is scalar window
% outputs:
%   bincent - centers of time bins
%   binvals - normalized likelihood (probability of z per y) for each bin

if nargin < 3 || isempty(trng)
    trng = max([abs(z(:)'); abs(y(:)')]); 
end

% ensure column vectors
y = y(:);
z = z(:);

% define bins: 1-ms resolution if times are in seconds, otherwise choose sensible bin
% choose bin width as median diff of combined times or default small value
allTimes = sort([y; z]);
if numel(allTimes) > 1
    dt = median(diff(allTimes));
    if dt <= 0
        dt = (2*trng)/100; % fallback
    end
else
    dt = (2*trng)/100;
end

% create bin edges covering [-trng, trng]
binedges = (dt/2):dt:trng; 
binedges = [-fliplr(binedges), binedges]; % force center at 0
bincent = (binedges(1:end-1) + binedges(2:end))/2; % should icl 0
nb = numel(bincent);
binvals = zeros(1, nb);

% for each y, compute relative times of z and histogram
% vectorized approach: for each y, consider z within window
for i = 1:numel(y)
    rel = z - y(i);
    sel = (rel >= -trng) & (rel <= trng);
    if any(sel)
        h = histcounts(rel(sel), binedges);
        binvals = binvals + h;
    end
end

% Normalize: convert counts to probability per y per bin
if isempty(y)
    binvals = zeros(size(binvals));
else
    binvals = binvals / numel(y) / dt; % rate (events per unit time) conditional on y
end
end