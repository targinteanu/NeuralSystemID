function Yp = myPredict(sys, tbl, k, showprog, usebuiltin)
% should work like predict with inputs system <sys>, input-output timetable
% data <tbl>, and prediction horizon <k> steps ahead 
% will call the subfunction that works for type of input <sys> 

% handle inputs 
if nargin < 5
    usebuiltin = false;
end
if nargin < 4
    showprog = true;
end
if nargin < 3
    k = 1;
end

if usebuiltin
    Yp = predictWrapper(sys, tbl, k);
else

% refer to best predict function 
c = class(sys);
if strcmpi(c, 'idss') || strcmpi(c, 'ss')
    Yp = myPredict2(sys, tbl, k, showprog);
elseif strcmpi(c, 'idNeuralStateSpace')
    Yp = myPredict1(sys, tbl, k, showprog);
elseif strcmpi(c, 'idpoly')
    Yp = myPredictAR(sys, tbl, k, showprog);
else
    warning(['System type ',c,' has not yet been tried.'])
    Yp = predictWrapper(sys, tbl, k);
end

end

    function Yp = predictWrapper(sys, tbl, k)
        k = k*tbl.Properties.TimeStep; % convert to seconds
        Yp = predict(sys, tbl, k);
        Yp.Time = Yp.Time + tbl.Time(1);
    end

end