function [hPlt, hAx] = myStackedPlot(tbl, vars)
% can be used like stackedplot() or plot() with one table <tbl> and
% variable(s) <vars>, but variables can be input as numeric indexes OR
% names, and names can be partial (i.e. 'AINP1' for 'AINP1   '). If no or
% empty variables, all columns of table will be plotted. 

evt = tbl.Properties.Events;

if nargin < 2
    vars = [];
end
if isempty(vars)
    vars = tbl.Properties.VariableNames;
end
ic = iscell(vars);

N = length(vars);
for c = 1:N
    if ic
        v = vars{c};
    else
        v = vars(c);
    end
    if isnumeric(v)
        v = tbl.Properties.VariableNames{v};
    else
        if ~sum(strcmp(tbl.Properties.VariableNames, v))
            iv = find(contains(tbl.Properties.VariableNames, v));
            if isempty(iv)
                error(['No table variables closely match ',v]);
            end
            if length(iv) > 1
                warning(['Multiple table variables closely match ',v,...
                    '; choosing ',tbl.Properties.VariableNames(iv(1))]);
            end
            iv = iv(1);
            v = tbl.Properties.VariableNames(iv);
        end
    end

    hAx(c) = subplot(N,1,c); 
    hPlt(c) = plot(tbl, v);
    grid on;
end

linkaxes(hAx, 'x');