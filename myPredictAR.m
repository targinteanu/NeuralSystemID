function Yp = myPredictAR(sys, tbl, k, showprog)
% should work like predict with inputs system <sys>, input-output timetable
% data <tbl>, and prediction horizon <k> steps ahead 
% works when <sys> is autonomous idpoly 

% setup 
if nargin < 4
    showprog = true;
end
progtick = .05; prog = 0; proggoal = .05;
if nargin < 3
    k = 1;
end
%k = k+1;
yp = nan(height(tbl), height(sys));
if iscell(sys.A)
    N = diag(cellfun(@width, sys.A)) - 1;
    N = unique(N);
    if length(N) > 1
        error('Channel AR models should be same order.')
    end
else
    N = width(sys.A) - 1;
end

% is autonomous? 
isAuton = ~width(sys); 
if ~isAuton
    error('Autonomous systems only.')
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

% first k+N predictions 
tblPad = array2timetable(zeros(N, width(tbl)), ...
    "SampleRate",tbl.Properties.SampleRate, ...
    "VariableNames",tbl.Properties.VariableNames);
tblPad.Time = tblPad.Time - tblPad.Time(end) + tbl.Time(1) - tbl.Properties.TimeStep;
TBL = [tblPad; tbl];
y1 = minisim(sys, TBL, 1, 2*N, k);
yp(1:(k+N),:) = y1((N+1):end,:);

% remaining predictions 
for ik2 = (N+1):(height(yp)-k)
    ik1 = ik2-N+1; ik3 = ik2+k;
    yy = minisim(sys, tbl, ik1, ik2, k);
    yp(ik3,:) = yy(end,:);
    prog = ik3/height(yp);
    if showprog
        if prog >= proggoal
            disp(['Simulating: ',num2str(100*prog),'% Complete'])
            proggoal = proggoal + progtick;
        end
    end
end

% convert output to timetable 
Yp = array2timetable(yp, "RowTimes",tbl.Time, ...
    "VariableNames",tbl.Properties.VariableNames);
Yp.Properties.VariableUnits = tbl.Properties.VariableUnits; 
% also copy events, userdata, customproperties, etc???
if ~isAuton
    Yp = [U, Yp]; 
end

% helper 
function ys = minisim(S, VT, i1, i2, k)
    % unsure if the non-autonomous case will work; do OutputTimes /
    % InitialCondition need to be duration/datetime/timetable to match U?
    vt = VT{i1:i2,:};
    ys = myFastForecastAR(S, vt, k);
    ys = [vt; ys];
end

end