function Yp = myPredict2(sys, tbl, k, showprog, fullstateavail)
% should work like predict with inputs system <sys>, input-output timetable
% data <tbl>, and prediction horizon <k> steps ahead 
% works when <sys> is idss

% setup 
if nargin < 5
    fullstateavail = true;
end
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

% estimate state 
if ~fullstateavail
    Y = tbl.Variables'; % row vectors 
    % Y = C * X
    C = sys.C;
    if (width(C) == height(C)) && (det(C) > 10*eps)
        Ci = C^-1;
    else
        %CTC = C' * C;
        %if det(CTC) > 10*eps
        %    Ci = (CTC ^ -1) * C';
        %else
            Ci = pinv(C);
        %end
    end
    % X ~~ Ci * Y
    X = Ci*Y;
    tbl = array2timetable(X', "RowTimes",tbl.Time);
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
    try
        ys = sim(S, UT, simopt);
    catch ME
        ys = sim(d2c(S), UT, simopt);
    end
end

end