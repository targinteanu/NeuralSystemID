function EvOut = AlphaOmegaTable2Depth(Tbl)
% 
% Input a table of alpha omega data with the Events property containing the
% file names. 
% Output a table with the depth, side (r/l), number, and file extracted
% from the file names. 
% 

Ev = Tbl.Properties.Events; 
lbl = Ev.EventLabels;

SIDE = repmat("", size(lbl));
N = nan(size(lbl));
DEPTH = nan(size(lbl));
FILE = nan(size(lbl));

for l = 1:height(lbl)
    [data, flag] = sscanf(lbl(l), '%ct%fd%ff%f %s');
    if flag > 4
        SIDE(l) = string(char(data(1)));
        N(l) = data(2);
        DEPTH(l) = data(3);
        FILE(l) = data(4);
    end
end

EvOut = timetable(Ev.Time, DEPTH, SIDE, N, FILE);
EvOut.Properties.VariableNames = {'DEPTH', 'SIDE', 'N', 'FILE'};

end