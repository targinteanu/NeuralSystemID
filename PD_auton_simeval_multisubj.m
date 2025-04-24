%% Parkinson's Disease (PD) project - evaluate multiple subjects' simulations
% Compare models with actual data across multiple subjects. Return plot of
% error and correlation over each prediction horizon. 
% Autonomous: no brain stimulation. 
% Only looking at AR and LTI models. 

%% user selects folder; necessary files are pulled 
folder = uigetdir; 
files = dir([folder,filesep,'*.mat']);

%% loop through each file 

CORS = {}; ERRS = {}; PCOR = {};
sysName = {'AR'; 'LTI'};
sysColr = {'b';  'r'};
fs = nan;

for f = files'
    load(fullfile(f.folder, f.name), 'cors', 'pcor', 'errs', 'dataTrain');

    fs_ = dataTrain.Properties.SampleRate; % Hz 
    if isnan(fs)
        fs = fs_;
    else
        if abs(fs - fs_) > 1e-4
            error('Incompatible sample rates.')
        end
    end

    % eliminate aggregate "channel" data 
    cors = cellfun(@(X) X(:,1:(end-1)), cors, 'UniformOutput',false)';
    pcor = cellfun(@(X) X(:,1:(end-1)), pcor, 'UniformOutput',false)';
    errs = cellfun(@(X) X(:,1:(end-1)), errs, 'UniformOutput',false)';

    CORS = [CORS, cors];
    ERRS = [ERRS, errs];
    PCOR = [PCOR, pcor];

    clear cors errs pcor dataTrain fs_
end

%% aggregate channels horizontally

cors = cell(height(CORS),1); 
errs = cors; pcor = cors;

% consistent length of horizons 
H = cellfun(@height, CORS); 
H = min(H,[],2);
for sys = 1:height(CORS)
    for subj = 1:width(CORS)
        cors{sys} = [cors{sys}, CORS{sys,subj}(1:H(sys),:)];
        errs{sys} = [errs{sys}, ERRS{sys,subj}(1:H(sys),:)];
        pcor{sys} = [pcor{sys}, PCOR{sys,subj}(1:H(sys),:)];
    end
end

%% stats and plotting 

lgd = sysName';
lgd = [lgd; repmat({'99% C.I.'}, size(lgd))];
lgd = lgd(:);

figure;

% corr 
ax(1) = subplot(2,1,1);
for sys = 1:height(cors)
    t = 1:H(sys); t = t-1; t = (t/fs); t = t';
    X = cors{sys};
    xMean = mean(X,2);
    xSEM = std(X,[],2) / sqrt(width(X));
    Tval = tinv(.005, width(X)-1);
    xCI = Tval*xSEM;
    plot(t, xMean, sysColr{sys}, 'LineWidth',2);
    hold on; grid on;
    patch([t; flipud(t)], [xMean+xCI; flipud(xMean-xCI)], ...
        sysColr{sys}, 'FaceAlpha',.25, 'EdgeColor','none');
end
%set(gca,'XScale','log');
ylabel('Pearsons \rho'); 
xlabel('Prediction Horizon (sec)');
legend(lgd, 'Location','northeast');

% err 
ax(2) = subplot(2,1,2);
for sys = 1:height(errs)
    t = 1:H(sys); t = t-1; t = (t/fs); t = t';
    X = errs{sys}*100;
    xMean = mean(X,2);
    xSEM = std(X,[],2) / sqrt(width(X));
    Tval = tinv(.005, width(X)-1);
    xCI = Tval*xSEM;
    plot(t, xMean, sysColr{sys}, 'LineWidth',2);
    hold on; grid on;
    patch([t; flipud(t)], [xMean+xCI; flipud(xMean-xCI)], ...
        sysColr{sys}, 'FaceAlpha',.25, 'EdgeColor','none');
end
%set(gca,'XScale','log');
ylabel('% RMSE'); 
xlabel('Prediction Horizon (sec)');
legend(lgd, 'Location','southeast');

linkaxes(ax, 'x');