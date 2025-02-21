function tbl = downsampleTimetable(TBL, rat) 
% Downsample an input <TBL> with sampling ratio <rat> such that the output
% table <tbl> will have one (averaged) entry for every <rat> input entries;
% <rat> should be an integer.  
if numel(TBL)
rat = round(rat); 
M = movmean(TBL.Variables, rat); 
m = downsample(M, rat);
t = downsample(TBL.Time, rat); 
tbl = array2timetable(m, "RowTimes",t, ...
    "VariableNames",TBL.Properties.VariableNames);
tbl.Properties.CustomProperties = TBL.Properties.CustomProperties; 
tbl.Properties.UserData = TBL.Properties.UserData; 
tbl.Properties.VariableUnits = TBL.Properties.VariableUnits;
tbl.Properties.VariableDescriptions = TBL.Properties.VariableDescriptions; 
tbl.Properties.Events = TBL.Properties.Events; % need to downsample this as well??
tbl.Properties.Description = [TBL.Properties.Description,' downsampled'];
else
    warning('No resampling done because table is empty!')
    tbl = TBL;
end