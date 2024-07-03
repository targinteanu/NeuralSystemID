function TT = ns2timetable(NS)

t = linspace(0,NS.MetaTags.DataPointsSec,NS.MetaTags.DataPoints);
t = seconds(t);
X = NS.Data'; X = single(X);
lbl = {NS.ElectrodesInfo.Label};
lbl = upper(lbl);

TT = array2timetable(X,"RowTimes",t,"VariableNames",lbl); 
%TT.Properties.SampleRate = NS.MetaTags.SamplingFreq;

% sample rate check 
fs_Signal = NS.MetaTags.SamplingFreq;
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