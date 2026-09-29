%% Artifact Removal Eval
% This code is designed to remove all stimulus artifacts from data that has
% been segmented. Artifact removal is causal in order to evaluate/compare
% techniques. 

%% obtain segmented data 

% file selection
thisfilename = mfilename("fullpath");
[fn,fp] = uigetfile('*SegmentData*.mat', 'Choose Segmented Data File');
SegmentedDataFullfile = fullfile(fp,fn);
load(SegmentedDataFullfile);
[~,fn,fe] = fileparts(fn);

% display details 
disp("Stimulus on channel(s): "); disp(channelNameStim);
disp("Trigger on channel(s): "); disp(channelNameTrig);
disp("Recording channel(s): "); disp(channelNameRec);

tblBL = tblsBaseline{1}; tblSt = tblsTrig{1};
fs = SampleRates(1);

stimtimes = tblSt.Properties.Events;
stimtimes = stimtimes(contains(stimtimes.EventLabels, 'Trig'), :);
stimtimes = stimtimes.Time;

%% (1) single-channel, single-band AR 

% define params
fbnd = [13,30]; % hz
ARord = 50;
trnLen = 500;
artDur = 16; % samples 
artStart = 1; % samples before
LMSstepsize = 0.5; 

% select chan
chname = channelNameRec(1); chidx = strcmp(tblBL.Properties.VariableNames, chname);
xBL = tblBL.(chname); xSt = tblSt.(chname);

% filter
BPF = buildFIRBPF(fs,fbnd(1),fbnd(2));
tblBLf = FilterTimetable(@(d,x) filtfilt(d,1,x), BPF, tblBL(:,chidx));
tblBLf = tblBLf(((-ceil(trnLen/2)+1):floor(trnLen/2)) + round(height(tblBL)/2), :);

% fit baseline model 
Mdl = ar(tblBLf,ARord,'yw'); 
wAR = -Mdl.A(2:end)/Mdl.A(1); wAR = fliplr(wAR); 

% use training residual to estimate process noise Q
xBLf = tblBLf{:,1}; 
xBLf_pred = myFastForecastAR(Mdl, xBLf(1:ARord), length(xBLf)-ARord);
xBLf = xBLf(ARord+1:end);
q = var(xBLf_pred - xBLf);

% setup LMS 
wLMS = zeros(artDur,1);
startbuffer = max(artDur, ARord);

% main loop
for stimtime = stimtimes
    [~,stimind] = min(abs(seconds(tblSt.Time - stimtime)));
    stimind = stimind - artStart;
    if stimind > startbuffer
    end
end