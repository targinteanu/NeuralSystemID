%% load data 
datasource = "/Users/torenarginteanu/Desktop/Data_PD/PD26N002/Neuro Omega/SPK_RT_SelectedTimes.mat"; 
waveclusfilename = "waveclusdata";
wcspikefilename = waveclusfilename+"_spikes.mat";
datasourcefolder = fileparts(datasource);

load(datasource)
data = TblD.CSPK_01;
spike_time = seconds(TblD.Time);
if any(diff(spike_time) < 0)
    error('time must be ascending and uniform.')
end
sr = 1/median(diff(spike_time)); % hz

tsel = (spike_time<(8440));
data = data(tsel);
spike_time = spike_time(tsel);

datadest = fullfile(datasourcefolder, waveclusfilename+".mat");
save(datadest, 'data', 'sr', '-v7.3');
Get_spikes(char(waveclusfilename+".mat"));
Do_clustering(char(wcspikefilename));

%% WaveClus spike detection directly from a timetable
addpath(genpath('/Users/torenarginteanu/Documents/MATLAB/wave_clus'));

% ---------------------------------------------------------
% WaveClus parameters
% ---------------------------------------------------------
par = set_parameters();

par.sr = sr;

% Detection band
par.detect_fmin = 300;
par.detect_fmax = 3000;
par.detect_order = 4;

% Threshold
par.stdmin = 5;
par.stdmax = 50;

% Spike polarity
par.detection = 'neg';      % 'pos', 'neg', or 'both'

% Samples retained around each detected spike
par.w_pre  = 20;
par.w_post = 44;

% Minimum time between detected spikes
par.ref_ms = 1.5;

% ---------------------------------------------------------
% Detect spikes
% ---------------------------------------------------------
[spikes, threshold, spike_index] = amp_detect(data(:)', par);

% spike_index is in samples
spike_time_sec = (spike_index - 1) / sr;

fprintf('Detected %d spikes.\n',length(spike_index));

% Put results into a convenient structure
spikeResults.spikes = spikes;
spikeResults.index = spike_index;
spikeResults.time = spike_time;
spikeResults.time_sec = spike_time_sec;
spikeResults.threshold = threshold;
spikeResults.par = par;