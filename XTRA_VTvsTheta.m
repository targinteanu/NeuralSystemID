filename = '/Users/torenarginteanu/Documents/Anderson Lab/XTRA/XTRA_neurophys-VT.xlsx'; % spreadsheet 
opts = detectImportOptions(filename);
opts.Sheet = 1;
opts.DataRange = 'A1';
data = readcell(filename, opts);
opts.Sheet = 2;
%opts = setvartype(opts, 'char'); % ensure flexible import
opts = setvartype(opts, 'double');
opts.VariableNamesRange = 1;      % first row contains variable names
opts.DataRange = 'A2';               % data starts on second row
anat = readtable(filename, opts);
anat = anat(:,2:6);

subjnames = string(anat.Properties.VariableNames);
subjrename = "Subject "+string(1:length(subjnames));

% reshape data into [<region> x <subj>] x <VT/BL/NB/BL2>
newsubjcol = cellfun(@ischar, data(1,:));
newsubjcol = newsubjcol | cellfun(@isstring, data(1,:));
for c = 1:length(newsubjcol)
    if newsubjcol(c)
        if ~contains(data{1,c}, 'XPD')
            newsubjcol(c) = false;
        end
    end
end
data2 = cell(length(subjnames), 4);
for s = 1:length(subjnames)
    subj = subjnames(s);
    c = find(strcmp(data(1,:), subj));
    if length(c) > 1
        error("Subject ID "+subj+" is not unique.");
    end
    if ~isempty(c)
        if c < length(newsubjcol)
            d = find(newsubjcol((c+1):end));
            if isempty(d)
                d = length(newsubjcol);
            else
                d = d(1)+c-1;
            end
        else
            d = length(newsubjcol);
        end
        for v = c:d
            vname = data{2,v};
            if strcmp(vname, 'VT')
                vi = 1;
            elseif strcmp(vname, 'BL')
                vi = 2;
            elseif strcmp(vname, 'BL2')
                vi = 3;
            elseif strcmp(vname, 'NB')
                vi = 4;
            else
                warning(['Unrecognized variable name ',vname]);
            end
            data2{s, vi} = [data{3:end, v}]';
        end
    end
end

figure; 
alignR = true;
VTcutoff = 15;
anatcutoff = 0.5;
mkr = {'s', 'd', 'x'};
clr = {[0.0660    0.4430    0.7450], ... blue 
       [0.8660    0.3290         0], ... red 
       [0.2310    0.6660    0.1960], ... green
       [0.5210    0.0860    0.8190], ... purple
       [0.4645    0.3470    0.0625], ... brown
       ...[0.9290    0.6940    0.1250], ... yellow
       [0.1840    0.7450    0.9370], ... teal
       [0.8190    0.0150    0.5450]  ... dark red
       };
for s = 1:length(subjnames)
    subj = subjnames(s);
    VT = data2{s,1};
    VTsel = VT > VTcutoff;
    anatsel = anat{:,s} > anatcutoff;
    if ~isempty(VT)
        for v = 2:4
            V = data2{s,v};
            if ~isempty(V)
                plot(VT, V, '.', ...
                    'Marker', mkr{v-1}, 'Color', clr{s}); 
                hold on;
                if sum(anatsel)
                    plot(VT(anatsel), V(anatsel), 'o', ...
                        'MarkerSize',10, 'Color',clr{s});
                end
                if sum(VTsel)
                    
Vsel_vals = V(VTsel);
VTsel_vals = VT(VTsel);
[p, fiteval] = polyfit(VTsel_vals, Vsel_vals, 1);
xfit = [VTcutoff, max(VTsel_vals)];
yfit = polyval(p, xfit);
plot(xfit, yfit, ':', 'Color', clr{s}, 'LineWidth', 1.5);
if alignR
    text(xfit(end), yfit(end), ['  y = ',num2str(p(1)),'x + ',num2str(p(2))], ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
        'Color',clr{s});
    text(xfit(end), yfit(end), ['  R^2 = ',num2str(fiteval.rsquared)], ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        'Color',clr{s});
else
    text(xfit(1), yfit(1), ['  y = ',num2str(p(1)),'x + ',num2str(p(2))], ...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom', ...
        'Color',clr{s});
    text(xfit(1), yfit(1), ['  R^2 = ',num2str(fiteval.rsquared)], ...
        'HorizontalAlignment','right', 'VerticalAlignment','top', ...
        'Color',clr{s});
end
alignR = ~alignR;

                end
            end
        end
    end
end
grid on;
xlabel('x = XTRA VT'); ylabel('y = Theta Power (dB)');