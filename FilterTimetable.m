function TTfilt = FilterTimetable(FiltFun, FiltObj, TTunfilt)
%
% Filter time-series data in a timetable 
%
% Inputs: 
%   FiltFun: function [e.g. @(d,x) filtfilt(d,x)] of form 
%            FilteredSignal = FiltFun(FiltObj, UnfilteredSignal) 
%            where FilteredSignal and UnfilteredSignal are matrix/vector 
%            columns over time.
%   FiltObj: digitalFilter object 
%   TTunfilt: unfiltered data in timetable format 
%
% Outputs: 
%   TTfilt: filtered data in timetable format 
%

%% sample rate check
fs_Filter = FiltObj.SampleRate; 
fs_Signal = TTunfilt.Properties.SampleRate;
Time_Signal = TTunfilt.Time; 
dT = seconds(diff(Time_Signal));
dTmax = max(dT); dTmin = min(dT);
fsmax = 1/dTmin; fsmin = 1/dTmax;
errthresh = .01;
if abs((fsmax-fs_Signal)/fs_Signal) > errthresh
    error('Timetable Sample Rate is incorrectly reported.')
end
if abs((fsmin-fs_Signal)/fs_Signal) > errthresh
    error('Timetable Sample Rate is incorrectly reported.')
end
if abs((fs_Filter-fs_Signal)/fs_Signal) > errthresh
    error('Filter and input signal have different sampling rates.')
end

%% filtering 
Xunfilt = table2array(TTunfilt);
Xfilt = FiltFun(FiltObj, Xunfilt);
TTfilt = array2timetable(Xfilt,"RowTimes",Time_Signal,...
    "VariableNames",TTunfilt.Properties.VariableNames); 
TTfilt.Properties.CustomProperties = TTunfilt.Properties.CustomProperties;

end