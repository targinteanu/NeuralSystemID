function Yp = myPredict2(sys, tbl, k, showprog)
% should work like predict with inputs system <sys>, input-output timetable
% data <tbl>, and prediction horizon <k> steps ahead 
% works when <sys> is idss

% setup 
if nargin < 4
    showprog = true;
end
progtick = .05; prog = 0; proggoal = .05;
if nargin < 3
    k = 1;
end
k = k+1;
yp = nan(height(tbl), height(sys));

% convert output to timetable 
Yp = array2timetable(yp, "RowTimes",tbl.Time, ...
    "VariableNames",tbl.Properties.VariableNames);
Yp.Properties.VariableUnits = tbl.Properties.VariableUnits; 
% also copy events, userdata, customproperties, etc???

% is autonomous? 
isAuton = ~width(sys); 
if ~isAuton
    tblInputInd = false(1, width(tbl));
    for NAME = sys.InputName
        name = NAME{:};
        tblInputInd = tblInputInd | strcmp(name, tbl.Properties.VariableNames);
    end
    %{
    tblInputInd = cellfun(...
        @(name) find(strcmp(name, tbl.Properties.VariableNames)), ...
        sys.InputName);
    %}
    U = tbl(:, tblInputInd);
    tbl = tbl(:, ~tblInputInd);
else
    U = [];
end

% first k predictions 
Yp(1:k,:) = minisim(sys, tbl, U, 1, k);

% remaining predictions 
for ik2 = (k+1):height(Yp)
    ik1 = ik2-k+1; 
    yy = minisim(sys, tbl, U, ik1, ik2);
    Yp(ik2,:) = yy(end,:);
    prog = ik2/height(Yp);
    if showprog
        if prog >= proggoal
            disp(['Simulating: ',num2str(100*prog),'% Complete'])
            proggoal = proggoal + progtick;
        end
    end
end

if ~isAuton
    Yp = [U, Yp]; 
end

% helper 
function ys = minisim(S, VT, UT, i1, i2)
    % unsure if the non-autonomous case will work; do OutputTimes /
    % InitialCondition need to be duration/datetime/timetable to match U?
    if ~isempty(UT)
        UT = UT(i1:i2,:);
    else
        UT = VT(i1:i2,[]);
    end
    simopt = simOptions(...
        'InitialCondition', VT{i1,:}');
    ys = sim(S, UT, simopt);
end

end