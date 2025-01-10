function [BandPower, SD, isOut, numOut, fig1] = tblChannelSummary(tbl, bandFreq)
% 
% Show some channel summary data, for example to inspect for bad channels. 
% 
% Inputs: 
%   tbl: timetable data (e.g. ECoG, EEG, etc) 
%   bandFreq: [lo, hi] freq of band of interest, e.g. [13, 30] for beta
%             band in Parkinson's Disease 
% 
% Outputs: 
%   BandPower: avg band power of each channel 
%   SD: standard deviation of each channel 
%   isOut: true/false matrix of outlier samples in each channel 
%   numOut: total # outliers of each channel 
%   fig1: figure handle

% TO DO:
% Should have a way to input outlier cutoffs and bypass outlier
% calculation. 
% might need to change this if multiple NS2 files with break in between,
% i.e. non-uniform sampling 

BandPower = bandpower(tbl.Variables, tbl.Properties.SampleRate, bandFreq);
SD = std(tbl);
isOut = isoutlier(tbl, 'mean');
numOut = sum(isOut);

fig1 = figure; 
subplot(4,1,1); stem(BandPower); axis tight;
ylabel('Band Power'); title('Channels Summary Data');
xticks(1:width(BandPower)); xticklabels(tbl.Properties.VariableNames);
subplot(4,1,2); stem(SD.Variables); axis tight;
ylabel('Standard Deviation'); 
xticks(1:width(SD)); xticklabels(tbl.Properties.VariableNames);
subplot(4,1,3); stem(numOut); axis tight;
ylabel('# Outliers'); 
xticks(1:width(numOut)); xticklabels(tbl.Properties.VariableNames);
subplot(4,1,4); boxplot(tbl.Variables, 'PlotStyle','compact', 'Symbol','.');
axis tight;
ylabel('Box Plot'); xlabel('Channel Name');
xticks(1:width(BandPower)); xticklabels(tbl.Properties.VariableNames);

end