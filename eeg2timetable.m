function TT = eeg2timetable(EEG)

t = milliseconds(EEG.times');
X = EEG.data'; 
lbl = {EEG.chanlocs.labels};
lbl = upper(lbl);

TT = array2timetable(X,"RowTimes",t,"VariableNames",lbl); 

end