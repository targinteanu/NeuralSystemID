function TT = eeg2timetable(EEG)

t = milliseconds(EEG.times');
X = EEG.data'; X = double(X);
lbl = {EEG.chanlocs.labels};
lbl = upper(lbl);

TT = array2timetable(X,"RowTimes",t,"VariableNames",lbl); 
%TT.Properties.SampleRate = EEG.srate;

% sample rate check 
fs_Signal = EEG.srate;
dT = seconds(diff(t));
dTmax = max(dT); dTmin = min(dT);
fsmax = 1/dTmin; fsmin = 1/dTmax;
errthresh = .01;
if abs((fsmax-fs_Signal)/fs_Signal) > errthresh
    warning('Sample Rate is incorrectly reported or inconsistent.')
end
if abs((fsmin-fs_Signal)/fs_Signal) > errthresh
    warning('Sample Rate is incorrectly reported or inconsistent.')
end

end