function T2 = myRemoveVars(T1, vars)

try
    T2 = removevars(T1, vars);
catch ME
    if contains(ME.identifier, 'UnrecognizedVarName')
        vars = string(vars);
        tblvars = T1.Properties.VariableNames;
        T2 = T1;
        for var = vars
            if sum(strcmp(tblvars, var))
                T2 = removevars(T2, var);
            else
                warning(['Table does not have variable ',char(var)]);
            end
        end
    else
        rethrow(ME)
    end
end

end