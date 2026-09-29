%% load data 
load("/Users/torenarginteanu/Desktop/Data_PD/PD26N002/Neuro Omega/SPK_RT_SelectedTimes.mat")
x = TblD.CSPK_01;
t = seconds(TblD.Time);
if any(diff(t) < 0)
    error('time must be ascending and uniform.')
end
fs = 1/median(diff(t)); % hz

tsel = (t<(8440));
x = x(tsel);
t = t(tsel);

addpath(genpath('/Users/torenarginteanu/Documents/MATLAB/wave_clus'));

% prepare data for wave_clus: expects a structure with .times or raw vector and sampling rate
% wave_clus main function is 'detect_spikes' in some versions; here use 'wavedet' and 'sortspikes' workflow
% Ensure signal is column vector
x = x(:)';

% create tmp file structure as expected by wave_clus spike detection functions
s_alg.sr = fs;
s_alg.datatype = 'continuous';
s_alg.spikes = [];
s_alg.detected = [];
% use default parameters from wave_clus: get default settings via 'get_default_settings' if available
if exist('get_default_settings','file')
    settings = get_default_settings;
else
    % minimal settings
    settings = struct();
    settings.stdmin = 5;    % threshold multiplier
    settings.w_pre = round(0.001*fs);  % 1 ms pre-peak
    settings.w_post = round(0.002*fs); % 2 ms post-peak
    settings.fmin = 300;
    settings.fmax = 3000;
end

% bandpass filter the signal as wave_clus expects (300-3000 Hz typical)
[b,a] = butter(3, [settings.fmin, settings.fmax]/(fs/2));
xf = filtfilt(b,a,double(x));

% detect spikes using wavedet if available
if exist('wavedet','file')
    % wavedet returns spikes and index of peaks
    [spikes, spikeTimes] = wavedet(xf, fs, settings);
    % ensure outputs consistent: spikes matrix columns are spike waveforms
    if size(spikes,1) > size(spikes,2)
        spikes = spikes';
    end
    s_alg.spikes = spikes;
    s_alg.detected = spikeTimes;
else
    % fallback simple threshold detection
    thr = settings.stdmin * median(abs(xf))/0.6745;
    above = xf > thr;
    d = diff([0 above 0]);
    starts = find(d==1);
    ends = find(d==-1)-1;
    peaks = zeros(1,length(starts));
    for k=1:length(starts)
        [~,idx] = max(xf(starts(k):ends(k)));
        peaks(k) = starts(k)+idx-1;
    end
    % extract waveforms
    w_pre = settings.w_pre;
    w_post = settings.w_post;
    L = w_pre + w_post + 1;
    spikes = zeros(L, numel(peaks));
    valid = peaks > w_pre & peaks <= (numel(xf)-w_post);
    peaks = peaks(valid);
    for k=1:numel(peaks)
        spikes(:,k) = xf(peaks(k)-w_pre:peaks(k)+w_post);
    end
    s_alg.spikes = spikes;
    s_alg.detected = peaks;
end

% optionally run sorting (if available)
if exist('do_clust','file') && ~isempty(s_alg.spikes)
    % do_clust expects spikes in columns
    try
        clu = do_clust(s_alg.spikes);
        s_alg.clu = clu;
    catch
        % ignore clustering errors
    end
end

% expose results
spike_waveforms = s_alg.spikes;
spike_sample_indices = s_alg.detected;

%% If wave_clus GUI functions are available, call them to display standard plots
if exist('gui_wave_clus','file') || exist('wave_clus','file') || exist('plot_all_waveforms','file') || exist('sortspikes','file')
    try
        % prepare temporary structure expected by wave_clus GUI/sorting functions
        wc_data = struct();
        % wave_clus versions differ; provide common fields:
        wc_data.spikes = s_alg.spikes;           % waveforms in columns
        wc_data.index = s_alg.detected;         % sample indices of spikes
        wc_data.sr = fs;                        % sampling rate
        % if clustering info exists, include
        if isfield(s_alg,'clu'); wc_data.clu = s_alg.clu; end
        % Try common GUI entry points in order of availability
        if exist('gui_wave_clus','file')
            % some versions accept a structure or workspace variables
            gui_wave_clus(wc_data);
        elseif exist('wave_clus','file')
            % older versions: wave_clus expects raw data file; try using sortspikes if available
            if exist('sortspikes','file')
                % call sortspikes with spikes and sampling rate if signature allows
                try
                    sortspikes(wc_data.spikes, wc_data.sr);
                catch
                    % fallback to plot_all_waveforms if available
                    if exist('plot_all_waveforms','file')
                        plot_all_waveforms(wc_data.spikes, []);
                    end
                end
            elseif exist('plot_all_waveforms','file')
                plot_all_waveforms(wc_data.spikes, []);
            end
        elseif exist('plot_all_waveforms','file')
            plot_all_waveforms(wc_data.spikes, []);
        else
            % last resort: if sortspikes exists, try calling with available args
            if exist('sortspikes','file')
                try
                    sortspikes(wc_data.spikes, wc_data.sr);
                catch
                    % silently continue if it fails
                end
            end
        end
    catch
        % if any GUI call fails, continue without throwing an error
    end
end
%% plot raw signal segment with detected spike times
figure('Name','Raw signal and detections','NumberTitle','off','Theme','light');
t_s = (0:numel(x)-1)/fs;
plot(t_s, x, 'k'); hold on;
if ~isempty(spike_sample_indices)
    plot(spike_sample_indices/fs, x(spike_sample_indices), 'ro', 'MarkerSize',6, 'LineWidth',1.2);
end
xlabel('Time (s)'); ylabel('Voltage'); title('Raw signal with detected spikes');
legend('Raw','Detections');

% plot average waveform and some individual waveforms
if ~isempty(spike_waveforms)
    figure('Name','Spike waveforms','NumberTitle','off','Theme','light');
    L = size(spike_waveforms,1);
    tw = (-settings.w_pre:settings.w_post)/fs*1000; % ms
    % plot subset of waveforms
    nplot = min(100, size(spike_waveforms,2));
    idx = round(linspace(1,size(spike_waveforms,2),nplot));
    plot(tw, spike_waveforms(:,idx), 'Color',[0.7 0.7 0.7]); hold on;
    plot(tw, mean(spike_waveforms,2), 'b', 'LineWidth',2);
    plot(tw, median(spike_waveforms,2), 'r--', 'LineWidth',1.5);
    xlabel('Time (ms)'); ylabel('Amplitude'); title('Extracted spike waveforms');
    legend('Individual','Mean','Median');
end

% if wave_clus clustering info available, show cluster results using built-ins
if isfield(s_alg,'clu') && ~isempty(s_alg.clu)
    % try to use plotAllWaveforms or similar if available
    if exist('plot_all_waveforms','file')
        figure('Name','Clustered waveforms (wave_clus)','NumberTitle','off','Theme','light');
        try
            plot_all_waveforms(s_alg.spikes, s_alg.clu);
        catch
            % fallback to simple colored scatter by cluster on PC space
            clu = s_alg.clu;
            clu = clu(:)';
            valid = clu>0 & clu<=max(clu);
            if any(valid)
                coeff = pca(double(s_alg.spikes)');
                scores = double(s_alg.spikes')*coeff(:,1:2);
                gscatter(scores(:,1), scores(:,2), clu, [], '.', 8);
                xlabel('PC1'); ylabel('PC2'); title('Spike clusters (PCA)');
            end
        end
    else
        % simple visualization: mean per cluster
        clu = s_alg.clu(:)';
        uniquec = unique(clu(clu>0));
        if ~isempty(uniquec)
            figure('Name','Cluster mean waveforms','NumberTitle','off','Theme','light'); hold on;
            cmap = lines(numel(uniquec));
            for i=1:numel(uniquec)
                c = uniquec(i);
                meanw = mean(s_alg.spikes(:,clu==c),2);
                plot(tw, meanw, 'Color',cmap(i,:), 'LineWidth',1.5);
            end
            xlabel('Time (ms)'); ylabel('Amplitude'); title('Mean waveform per cluster');
            legend(arrayfun(@(c) sprintf('Cluster %d',c), uniquec, 'UniformOutput',false));
        end
    end
end

% raster / PSTH style: plot spike times as raster and compute rate
if ~isempty(spike_sample_indices)
    figure('Name','Spike times and firing rate','NumberTitle','off','Theme','light');
    subplot(2,1,1);
    stem(spike_sample_indices/fs, ones(size(spike_sample_indices)), 'Marker','none');
    xlim([t_s(1) t_s(end)]);
    ylabel('Spikes'); title('Spike raster');
    subplot(2,1,2);
    binw = 0.1; % 100 ms
    edges = t_s(1):binw:t_s(end);
    counts = histcounts(spike_sample_indices/fs, edges);
    bar(edges(1:end-1)+binw/2, counts/binw, 'k');
    xlabel('Time (s)'); ylabel('Firing rate (Hz)'); title(sprintf('Firing rate (bin=%.0f ms)',binw*1000));
end