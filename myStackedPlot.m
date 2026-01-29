function [hPlt, hAx] = myStackedPlot(tbl, vars, IsOutlier)
% can be used like stackedplot() or plot() with one table <tbl> and
% variable(s) <vars>, but variables can be input as numeric indexes OR
% names, and names can be partial (i.e. 'AINP1' for 'AINP1   '). If no or
% empty variables, all columns of table will be plotted. 

evt = tbl.Properties.Events;

if nargin < 3
    IsOutlier = [];
end
if nargin < 2
    vars = [];
end
if isempty(vars)
    vars = tbl.Properties.VariableNames;
end
ic = iscell(vars);

N = length(vars);
for c = 1:N

    % what variable to plot
    if ic
        v = vars{c};
    else
        v = vars(c);
    end
    if isnumeric(v)
        iv = v;
        v = tbl.Properties.VariableNames{v};
    else
        v = char(v);
        if ~sum(strcmp(tbl.Properties.VariableNames, v))
            iv = find(contains(tbl.Properties.VariableNames, v));
            if isempty(iv)
                error(['No table variables closely match ',char(v)]);
            end
            if length(iv) > 1
                warning(['Multiple table variables closely match ',v,...
                    '; choosing ',tbl.Properties.VariableNames{iv(1)}]);
            end
            iv = iv(1);
            v = tbl.Properties.VariableNames(iv);
        end
        iv = find(strcmp(tbl.Properties.VariableNames, v));
        iv = iv(1); 
    end
        if isempty(tbl.Properties.VariableUnits)
            u = '';
        else
            u = tbl.Properties.VariableUnits{iv};
        end
        if isempty(tbl.Properties.VariableDescriptions)
            d = '';
        else
            d = tbl.Properties.VariableDescriptions{iv};
            if length(d) > 20 % adjust max length here
                d = [d(1:20),'...'];
            end
            di = find(d==' '); % comma instead?
            if ~isempty(di)
                d = [d(1:di),'...'];
            end
        end
        if ~isempty(IsOutlier)
            io = IsOutlier(:,iv);
        end
    hAx(c,1) = subplot(N,1,c); 
    x = tbl.Time;
    y = tbl.(v);
    yl = [min(y(:)), max(y(:))];

    % plot events
    if ~isempty(evt)
        hold on;
        %yl = ylim();
        xx = evt.Time; 
        yy1 = yl(2)*ones(size(xx));
        yy2 = yl(1)*ones(size(xx));
        hPlt(c,2) = stem(xx, yy1, 'Marker','none', 'Color',[.5,.5,.5], 'LineWidth',.1);
        hPlt(c,3) = stem(xx, yy2, 'Marker','none', 'Color',[.5,.5,.5], 'LineWidth',.1);
        %ylim(yl);
    end

    % plot timed data
    hold on;
    hPlt(c,1) = plot(x, y);

    % plot outliers 
    if ~isempty(IsOutlier)
        if sum(io)
            hold on;
            hPlt(c,4) = plot(x(io), y(io), 'o');
        end
    end

    % axis labeling 
    grid on; axis tight;
    if c < N
        xticklabels([]);
    else
        xlabel('Time');
    end
    ylabel({v,['(',u,')'],d});

end

linkaxes(hAx, 'x');