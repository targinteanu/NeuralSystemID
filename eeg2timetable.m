function TT = eeg2timetable(EEG)

t = milliseconds(EEG.times');
X = EEG.data'; 
lbl = {EEG.chanlocs.labels};

TT = array2timetable(X,"RowTimes",t,"VariableNames",lbl); 

end