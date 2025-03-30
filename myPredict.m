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

%isfield_ = @(sys, fld) sum(strcmp(fieldnames(sys), fld));
%if isfield_(sys, 'StateName') && isfield_(sys, 'OutputName')
try
fullStateAvail = isequal(size(sys.StateName), size(sys.OutputName));
    % TO DO: should be a more robust check if full state output; maybe if
    % names are the same (not just size)? 
%else
catch
    fullStateAvail = true;
end

if usebuiltin
    Yp = predictWrapper(sys, tbl, k);
else

% refer to best predict function 
c = class(sys);
if strcmpi(c, 'idss') || strcmpi(c, 'ss')
    if fullStateAvail
        %Yp = myPredict2(sys, tbl, k, showprog, fullStateAvail);
        if showprog || ~usebuiltin
            warning('myPredict has been programmed to use built-in for this type.')
        end
        Yp = predict(sys, tbl, k, ...
            predictOptions('InitialCondition',tbl{1,:}')); 
        Yp.Time = Yp.Time + tbl.Time(1);
    else
        if showprog || ~usebuiltin
            warning('myPredict has been programmed to use built-in for this type.')
        end
        Yp = predict(sys, tbl, k, ...
            predictOptions('InitialCondition','z')); 
        Yp.Time = Yp.Time + tbl.Time(1);
    end
elseif strcmpi(c, 'idnlhw')
    if showprog || ~usebuiltin
        warning('myPredict has been programmed to use built-in for this type.')
    end
    Yp = predict(sys, tbl, k, ...
        predictOptions('InitialCondition','z'));
    Yp.Time = Yp.Time + tbl.Time(1);
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