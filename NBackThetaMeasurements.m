%% Calculate and display Theta measurements during N-Back task 

% Loading ECoG data into timetables sampled at 1kHz. 
%   BL = baseline data; NB = N-Back cognitive task 
bl = openNSx('PD25N009/baseline001.ns2', 'uV'); BL = ns2timetable(bl);
nb = openNSx('PD25N009/nback002.ns2', 'uV'); NB = ns2timetable(nb);

%% Looping through the channels 
% 
% Go through each channel. Filter into the theta and gamma ranges using
% filtfilt. Use calcPAC to get the theta-gamma coupling. 
% 

SamplingFreq = ns2.MetaTags.SamplingFreq;
t = NS2.Time';
tRel = seconds(t-t(1));
channelNames = NS2.Properties.VariableNames;

% setup bandpass filters
bpfTheta = buildFIRBPF(SamplingFreq, 4, 9); % Theta: 4 to 9 Hz
bpfGamma = buildFIRBPF(SamplingFreq, 30, 100); % low Gamma 

for channelIndex = 1:width(BL)

channelName = BL.Properties.VariableNames{channelIndex}
dataOneChannel = BL{:,channelIndex}';

% filtfilt ...
% calcPAC ...

end

%% show bar plot of each channel's PAC in baseline vs nback

figure;
% bar ...