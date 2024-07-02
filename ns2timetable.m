function TT = ns2timetable(NS)

t = linspace(0,NS.MetaTags.DataPointsSec,NS.MetaTags.DataPoints);
t = seconds(t);
X = NS.Data'; X = single(X);
lbl = {NS.ElectrodesInfo.Label};
lbl = upper(lbl);

TT = array2timetable(X,"RowTimes",t,"VariableNames",lbl); 

end