filename = '/Users/torenarginteanu/Documents/Anderson Lab/XTRA/XTRA_neurophys-VT.xlsx'; % spreadsheet 
%{
opts = detectImportOptions(filename);
opts.Sheet = 3;
opts.DataRange = 'A1';
%}
theta = readcell(filename, 'Sheet',3, 'Range','A:S');

opts = detectImportOptions(filename);
opts.Sheet = 1;
opts = setvartype(opts, 'double');
opts.VariableNamesRange = 1;      % first row contains variable names
opts.DataRange = 'A2';               % data starts on second row
anatAH = readtable(filename, opts);
anatAH = anatAH(:,2:6);

opts = detectImportOptions(filename);
opts.Sheet = 2;
opts = setvartype(opts, 'double');
opts.VariableNamesRange = 1;      % first row contains variable names
opts.DataRange = 'A2';               % data starts on second row
anatHCP = readtable(filename, opts);
anatHCP = anatHCP(:,2:6);

subjnames = string(anatHCP.Properties.VariableNames);
subjrename = "Subject "+string(1:length(subjnames));

%% reshape data into [<region> x <subj>] x <VT/BL/NB/BL2>
newsubjcol = cellfun(@ischar, theta(1,:));
newsubjcol = newsubjcol | cellfun(@isstring, theta(1,:));
for c = 1:length(newsubjcol)
    if newsubjcol(c)
        if ~contains(theta{1,c}, 'XPD')
            newsubjcol(c) = false;
        end
    end
end
data = cell(length(subjnames), 4);
for s = 1:length(subjnames)
    subj = subjnames(s);
    c = find(strcmp(theta(1,:), subj));
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
            vname = theta{2,v};
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
            data{s, vi} = [theta{3:end, v}]';
        end
    end
end

%%
clr = {[0.0660    0.4430    0.7450], ... blue 
       [0.8660    0.3290         0], ... red 
       [0.2310    0.6660    0.1960], ... green
       [0.5210    0.0860    0.8190], ... purple
       [0.4645    0.3470    0.0625], ... brown
       ...[0.9290    0.6940    0.1250], ... yellow
       [0.1840    0.7450    0.9370], ... teal
       [0.8190    0.0150    0.5450]  ... dark red
       };
figure; 
for s = 1:length(subjnames)
    subj = subjnames(s);
    VT = data{s,1};
    aHCP = anatHCP{:,s}; 
    aAH = anatAH{:,s};
    if ~isempty(VT)
        %plot3(aHCP, aAH, VT, 's', 'Color',clr{s}, 'LineWidth',1.5);
        plot(.5*(aHCP+aAH), VT, 's', 'Color',clr{s}, 'LineWidth',1.5);
        hold on;
    end
end
grid on;
%xlabel('HCP'); ylabel('Al-Hakim'); zlabel('XTRA VT');
xlabel('anatomy'); ylabel('XTRA VT');
title('Anatomy-VT Comparison');
legend(anatHCP.Properties.VariableNames, 'Location','eastoutside');

%%
figure; 
alignR = true;
%VTcutoff = 15;
anatcutoff = 0.01;
mkr = {'s', 'o', 'x'};
for s = 1:length(subjnames)
    subj = subjnames(s);
    VT = data{s,1};
    %VTsel = VT > VTcutoff;
    anatsel = anatHCP{:,s} > anatcutoff;
    if ~isempty(VT)
        VT = VT(anatsel);
        for v = 2:4
            V = data{s,v};
            if ~isempty(V)
                V = V(anatsel);
                plot(VT, V, '.', ...
                    'Marker', mkr{v-1}, 'Color', clr{s}, 'LineWidth',1.5); 
                hold on;
                %if sum(VTsel)
                    
%Vsel_vals = V(VTsel);
%VTsel_vals = VT(VTsel);
Vsel_vals = V; VTsel_vals = VT;
[p, fiteval] = polyfit(VTsel_vals, Vsel_vals, 1);
rho = corr(VTsel_vals, Vsel_vals, "Type","Spearman");
xfit = [min(VTsel_vals), max(VTsel_vals)];
yfit = polyval(p, xfit);
plot(xfit, yfit, ':', 'Color', clr{s}, 'LineWidth', 1.5);
txt1 = ['  \it y\rm = ',num2str(p(1),3),'\it x\rm + ',num2str(p(2),3),'  '];
%txt2 = ['  R^2 = ',num2str(fiteval.rsquared,3),'  '];
txt2 = ['  \rho = ',num2str(rho,3),'  '];
if alignR
    text(xfit(end), yfit(end), txt1, ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
        'Color',clr{s});
    text(xfit(end), yfit(end), txt2, ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        'Color',clr{s});
else
    text(xfit(1), yfit(1), txt1, ...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom', ...
        'Color',clr{s});
    text(xfit(1), yfit(1), txt2, ...
        'HorizontalAlignment','right', 'VerticalAlignment','top', ...
        'Color',clr{s});
end
alignR = ~alignR;

                %end
            end
        end
    end
end
grid on;
xlabel('\it x\rm = XTRA VT'); ylabel('\it y\rm = Theta Power (dB)');
xlim([5 25]);