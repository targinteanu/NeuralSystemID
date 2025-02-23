%% user selects folder; necessary files are pulled 
folder = uigetdir; 
files = dir([folder,filesep,'*.mat']);

%% loop through each file 

hznsAll = []; 
errsAll = {}; corsAll = {}; pcorsAll = {};
namesAll = {};
for f = files'
    load(fullfile(f.folder, f.name), 'hzns', 'errsTest', 'corsTest', 'pcorsTest', 'dataTest');
    % row = hzn 
    % col = channel 
    % sheet = model 

    errsAll = [errsAll, errsTest]; 
    corsAll = [corsAll, corsTest]; 
    pcorsAll = [pcorsAll, pcorsTest];
    
    chNames = dataTest.Properties.VariableNames;
    chNames = cellfun(@(n) n(1:min(4, length(n))), chNames, ... -> 'LSxx
        'UniformOutput',false);
    pName = f.name; 
    pName = pName(1:min(8, length(pName))); % -> PDyyNxxx
    namesAll = [namesAll, [pName; {chNames}]];

    fs = dataTest.Properties.SampleRate; % Hz
    if isnan(fs)
        fs = 1/seconds(mean(diff(dataTest.Time)));
    end
    hzns = round(hzns*1000/fs); % ms
    hznsAll = [hznsAll; hzns]; % horizontal elements 

    clear dataTest hzns errsTest corsTest pcorsTest
end

hznsAll = unique(hznsAll, "rows");
if height(hznsAll) > 1
    warning('Prediction horizons may be inconsistent between subjects.')
    hznsAll
    hznsAll = mode(hznsAll);
end

pcorsAllCat = []; corsAllCat = []; errsAllCat = [];
for subj = 1:width(pcorsAll)
    pcorsAllCat = [pcorsAllCat, pcorsAll{subj}];
    corsAllCat = [corsAllCat, corsAll{subj}];
    errsAllCat = [errsAllCat, errsAll{subj}];
end

%% num of stat sig correlated channels 
alph = .05; % confidence level, uncorrected 

nChanAll = width(pcorsAllCat);
alph_ = alph/nChanAll; % Bonferroni correction 
hAll = pcorsAllCat < alph_;
hCount = sum(hAll, 2);

for h = 1:height(hCount)
    disp(['At the ',num2str(hznsAll(h)),'-ms-ahead prediction, ' ...
        'the AR model was positively correlated at ',num2str(hCount(h,1,1)) ...
        ' and the LTI model was positively correlated at ',num2str(hCount(h,1,2)) ...
        ' out of ',num2str(nChanAll),' total channels (p < ', ...
        num2str(alph_,1),').'])
end

%% bar plots

errsMean = mean(100*errsAllCat, 2); errsSD = std(100*errsAllCat, [], 2);
corsMean = mean(corsAllCat, 2); corsSD = std(corsAllCat, [], 2);
errsP = ones(length(hznsAll), 1); corsP = errsP;
H = length(hznsAll);
for h = 1:H
    [~,errsP(h,1)] = ttest(errsAllCat(h,:,1), errsAllCat(h,:,2));
    [~,corsP(h,1)] = ttest(corsAllCat(h,:,1), corsAllCat(h,:,2));
end

figure; 
mkrs = {'o', 'x', 'd', '+', 's', 'v', '^'};
W = .6; % width range for individual points scatter 
for h = 1:H

    % RMSE 
    ax(h,1) = subplot(H,2, 2*(h-1)+1); 
    b = bar(squeeze(errsMean(h,:,:)), 'FaceAlpha',.2, 'LineWidth',2); grid on; hold on;
    x = b.XData; y = b.YData;  
    axis tight;

    % individual pts 
    for subj = 1:width(namesAll)
        Y = 100*errsAll{subj}(h,:,:);
        yy = []; xx = [];
        for m = 1:size(Y,3)
            yyy = Y(:,:,m); yy = [yy, yyy]; 
            xxx = x(m)*ones(size(yyy));
            xxx = xxx + ( (0:(length(xxx)-1))*W/length(xxx) -W/2 ); % distribute points horizontally
            xx = [xx, xxx];
        end
        plot(xx, yy, [mkrs{subj},'r'], 'LineWidth',1.5);
    end
    errorbar(x,y, squeeze(errsSD(h,:,:)),squeeze(errsSD(h,:,:)), '.k', 'LineWidth',4);

    % pval 
    p = errsP(h);
    if p < alph
        yl = ylim(); yl_ = .8*max(yl);
        errorbar(mean(x),yl_, 0,0, diff(x)/2,diff(x/2), '.k', 'LineWidth',2);
        text(mean(x),yl_, ['p = ',num2str(p,1)], ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize',13);
    end
    legend(['Mean', namesAll(1,:), '±1SD'])
    xticks(x); xticklabels({'AR', 'LTI'}); xlabel('Model');
    ylabel('% RMSE');
    title([num2str(hznsAll(h)),'ms prediction'])
    set(gca, 'FontSize',16)

    % corr 
    ax(h,1) = subplot(H,2, 2*(h-1)+2); 
    b = bar(squeeze(corsMean(h,:,:)), 'FaceAlpha',.2, 'LineWidth',2); grid on; hold on;
    x = b.XData; y = b.YData; 
    axis tight;

    % individual pts 
    for subj = 1:width(namesAll)
        Y = corsAll{subj}(h,:,:);
        yy = []; xx = [];
        for m = 1:size(Y,3)
            yyy = Y(:,:,m); yy = [yy, yyy]; 
            xxx = x(m)*ones(size(yyy));
            xxx = xxx + ( (0:(length(xxx)-1))*W/length(xxx) -W/2 ); % distribute points horizontally
            xx = [xx, xxx];
        end
        plot(xx, yy, [mkrs{subj},'r'], 'LineWidth',1.5);
    end
    errorbar(x,y, squeeze(corsSD(h,:,:)),squeeze(corsSD(h,:,:)), '.k', 'LineWidth',4); 

    % pval 
    p = corsP(h);
    if p < alph
        yl = ylim(); yl_ = .8*max(yl);
        errorbar(mean(x),yl_, 0,0, diff(x)/2,diff(x/2), '.k', 'LineWidth',2);
        text(mean(x),yl_, ['p = ',num2str(p,1)], ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize',13);
    end
    legend(['Mean', namesAll(1,:), '±1SD'])
    xticks(x); xticklabels({'AR', 'LTI'}); xlabel('Model');
    ylabel('Pearsons \rho');
    title([num2str(hznsAll(h)),'ms prediction'])
    set(gca, 'FontSize',16)

end