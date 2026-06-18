subjIDdict = {...
    'PD23N005', 'XPD052'; ...
    'PD23N009', 'XPD053'; ...
    'PD24N003', 'XPD056'; ...
    'PD24N010', 'XPD058'; ...
    'PD24N008', 'XPD059'; ...
    'PD25N001', 'XPD061'; ...
    'PD25N009', 'XPD062'; ...
    'PD25N008', 'XPD063'; ...
    'PD25N010', 'XPD064'; ...
    'PD25N013', 'XPD065'; ...
    'PD26N002', 'XPD066'};
subjIDdict = containers.Map(subjIDdict(:,1), subjIDdict(:,2));

filename = '/Users/torenarginteanu/Documents/Anderson Lab/XTRA/dataL_91f313d618c254/Rows_master.xlsx'; % spreadsheet 

% VT
opts = detectImportOptions(filename);
opts.Sheet = 1;
opts = setvartype(opts, 'double');
opts.VariableNamesRange = 1;      % first row contains XPD names
opts.DataRange = 'A3';               % data starts on 3rd row
tblVT = readtable(filename, opts);

% PAC, power 
tblEph = readtable(filename, "Sheet",2, "Range","A1:CU22", "ReadVariableNames",true);
%{
opts = detectImportOptions(filename);
opts.Sheet = 2;                   % low theta X low gamma
opts = setvartype(opts, 'double');
opts.VariableNamesRange = 1;      % first row contains description
opts.DataRange = 'A2:CU22';               % data starts on 2nd row
tblEph = readtable(filename, opts);
%}

filename = '/Users/torenarginteanu/Documents/Anderson Lab/XTRA/XTRA_neurophys-VT.xlsx'; % spreadsheet 

% Al-Hakim
opts = detectImportOptions(filename);
opts.Sheet = 1;
opts = setvartype(opts, 'double');
opts.VariableNamesRange = 1;      % first row contains variable names
opts.DataRange = 'A2';               % data starts on second row
anatAH = readtable(filename, opts);
anatAH = anatAH(:,2:6);

% HCP 
opts = detectImportOptions(filename);
opts.Sheet = 2;
opts = setvartype(opts, 'double');
opts.VariableNamesRange = 1;      % first row contains variable names
opts.DataRange = 'A2';               % data starts on second row
anatHCP = readtable(filename, opts);
anatHCP = anatHCP(:,2:6);

subjnames = string(anatHCP.Properties.VariableNames);
subjrename = "Subject "+string(1:length(subjnames));

%% parse ephys variable names 

varnames = tblEph.Properties.VariableNames;
tblEphProp = cell(4, length(varnames));
    % 1: PD#
    % 2: rec type 
    % 3: band (power) or PAC 
    % 4: avg/std/#

for v = 1:length(varnames)
    tblEphProp(:,v) = strsplit(varnames{v}, '_')';
end

for v = 1:length(varnames)

    % rename subj PD# -> XPD#
    tblEphProp{1,v} = subjIDdict(tblEphProp{1,v});

    % rename band/PAC
    bnd = tblEphProp{3,v};
    if strcmpi(bnd,'The')
        tblEphProp{3,v} = 'Theta Power (dB)';
    elseif strcmpi(bnd,'Gam')
        tblEphProp{3,v} = 'Gamma Power (dB)';
    end

    % rename rec type -> BL/NB/BL2
    rectype = tblEphProp{2,v};
    if strcmpi(rectype, 'Baseline')
        tblEphProp{2,v} = 'BL';
    elseif contains(lower(rectype), 'nback')
        tblEphProp{2,v} = 'NB';
    else 
        tblEphProp{2,v} = 'BL2';
    end
end

%% reshape data 
% into <subj>[<row> x <pwr/PAC> x <val/errorbar>] x <VT/BL/NB/BL2>
data = cell(length(subjnames), 4);

bandnames = unique(tblEphProp(3,:));
vnames =["VT", "VT"; 
         "BL", "Rest Baseline"; 
         "BL2", "Rest After Task"; 
         "NB", "Task"];
vnamesdict = containers.Map(vnames(:,1), vnames(:,2));

for s = 1:length(subjnames)
    subj = subjnames(s);
    c = find(strcmp(tblEphProp(1,:), subj));
    if ~isempty(c)
        tblEphProp_s = tblEphProp(:,c); tblEph_s = tblEph{:,c};

        % VT: v=1 
        data{s,1} = tblVT.(subj);

        for v = 2:height(vnames) % BL, BL2, NB
            vname = vnames(v,1);
            cc = find(strcmp(tblEphProp_s(2,:), vname));
            tblEphProp_sv = tblEphProp_s(:,cc); tblEph_sv = tblEph_s(:,cc);
            data_sv = nan(height(tblEph),length(bandnames),2);
            for b = 1:length(bandnames)
                bandname = bandnames{b};
                ccc = find(strcmp(tblEphProp_sv(3,:), bandname));
                tblEphProp_svb = tblEphProp_sv(:,ccc); tblEph_svb = tblEph_sv(:,ccc);
                cX = find(strcmp(tblEphProp_svb(4,:), 'AVG'));
                cV = find(strcmp(tblEphProp_svb(4,:), 'STD'));
                cN = find(strcmp(tblEphProp_svb(4,:), 'NUM'));
                X = tblEph_svb(:,cX); V = tblEph_svb(:,cV); N = tblEph_svb(:,cN);
                V = tinv(.975, N-1).*(V./sqrt(N)); % STD -> 95% CI
                if ~isempty(X)
                    data_sv(:,b,1) = X; % Store average values
                end
                if ~isempty(V)
                    data_sv(:,b,2) = V; % Store errorbar values
                end
            end
            data{s,v} = data_sv;
        end
    end
end

%%
clr = {[0.0660    0.4430    0.7450], ... blue 
       [0.8660    0.3290         0], ... red 
       [0.2310    0.6660    0.1960], ... green
       [0.5210    0.0860    0.8190], ... purple
       [0.6193    0.4627    0.0833], ... gold
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

%% Power vs VT: all in one 

for b = 1:length(bandnames)
bandname = bandnames{b};

figure; 
lgd = [];
alignR = true;
%VTcutoff = 15;
anatcutoff = 0.01;
mkr = {'s',  'o',  'x'; 
       '--', '-.', ':'};

for s = 1:length(subjnames)

    subj = subjnames(s);
    VT = data{s,1};
    %VTsel = VT > VTcutoff;
    anatsel = anatHCP{:,s} > anatcutoff;
    if ~isempty(VT)
        VT = VT(anatsel);
        for v = 2:height(vnames)
            g = 1 + 2*(v-2)/3; % gamma color brightness adjustment
            D = data{s,v};
            if ~isempty(D)
                D = D(anatsel,b,:);
                Derr = D(:,2); D = D(:,1);
                errorbar(VT, D, Derr, -Derr, '.', ...
                    'Marker', mkr{1,v-1}, 'Color', clr{s}.^g, 'LineWidth',1.5); 
                hold on;
                lgd = [lgd, subj+":", vnames(v,2)];
                %if sum(VTsel)
                    
%Vsel_vals = V(VTsel);
%VTsel_vals = VT(VTsel);
Vsel_vals = D; VTsel_vals = VT;
[p, fiteval] = polyfit(VTsel_vals, Vsel_vals, 1);
rho = corr(VTsel_vals, Vsel_vals, "Type","Spearman");
xfit = [min(VTsel_vals), max(VTsel_vals)];
yfit = polyval(p, xfit);
plot(xfit, yfit, mkr{2,v-1}, 'Color', clr{s}.^g, 'LineWidth', 1.5);
txt1 = ['  \it y\rm = ',num2str(p(1),'%.1f'),'\it x\rm + ',num2str(p(2),'%.1f'),'  '];
%txt2 = ['  R^2 = ',num2str(fiteval.rsquared,'%+.2f'),'  '];
txt2 = ['  \rho = ',num2str(rho,'%+.2f'),'  '];
if alignR
    text(xfit(end), yfit(end), txt1, ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
        'Color',clr{s}.^g);
    text(xfit(end), yfit(end), txt2, ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        'Color',clr{s}.^g);
else
    text(xfit(1), yfit(1), txt1, ...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom', ...
        'Color',clr{s}.^g);
    text(xfit(1), yfit(1), txt2, ...
        'HorizontalAlignment','right', 'VerticalAlignment','top', ...
        'Color',clr{s}.^g);
end
alignR = ~alignR;

                %end
            end
        end
    end
end

grid on;
xlabel('\it x\rm = XTRA VT'); ylabel(['\it y\rm = ',bandname]);
xlim([6 24]);
legend(lgd, 'Location','eastoutside')

end

%% Power vs VT: baseline all subjects 

figure; 
alignR = true;

% first: all scatter plot
for s = 1:length(subjnames)
    subj = subjnames(s);
    VT = data{s,1};
    anatsel = anatHCP{:,s} > anatcutoff;
    if ~isempty(VT)
        VT = VT(anatsel);
        v = 2; % baseline 
        D = data{s,v};
        if ~isempty(D)
            D = D(anatsel);
            plot(VT, D, '.', ...
                'Marker', mkr{1,v-1}, 'Color', clr{s}, 'LineWidth',1.5); 
            hold on;
        end
    end
end


% second: all trendlines
for s = 1:length(subjnames)
    subj = subjnames(s);
    VT = data{s,1};
    anatsel = anatHCP{:,s} > anatcutoff;
    if ~isempty(VT)
        VT = VT(anatsel);
        v = 2; % baseline 
        D = data{s,v};
        if ~isempty(D)
            D = D(anatsel);

Vsel_vals = D; VTsel_vals = VT;
[p, fiteval] = polyfit(VTsel_vals, Vsel_vals, 1);
rho = corr(VTsel_vals, Vsel_vals, "Type","Spearman");
xfit = [min(VTsel_vals), max(VTsel_vals)];
yfit = polyval(p, xfit);
plot(xfit, yfit, mkr{2,v-1}, 'Color', clr{s}, 'LineWidth', 1.5);
txt1 = ['  \it y\rm = ',num2str(p(1),'%.1f'),'\it x\rm + ',num2str(p(2),'%.1f'),'  '];
txt2 = ['  \rho = ',num2str(rho,'%+.2f'),'  '];
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

        end
    end
end


grid on;
xlabel('\it x\rm = XTRA VT'); ylabel(['\it y\rm = ',bandname]);
xlim([6 24]);
legend(anatHCP.Properties.VariableNames, 'Location','eastoutside')
title('Baseline All Subjects')

%% Power vs VT: all cond each subject

for s = 1:length(subjnames)
figure; 
lgd = [];
alignR = true;

    subj = subjnames(s);
    VT = data{s,1};
    anatsel = anatHCP{:,s} > anatcutoff;
    if ~isempty(VT)
        VT = VT(anatsel);
        for v = 2:4
            g = 1 + 2*(v-2)/3; % gamma color brightness adjustment
            D = data{s,v};
            if ~isempty(D)
                D = D(anatsel);
                plot(VT, D, '.', ...
                    'Marker', mkr{1,v-1}, 'Color', clr{s}.^g, 'LineWidth',1.5); 
                hold on;
                lgd = [lgd, subj+":", vnames(v,2)];
                    
Vsel_vals = D; VTsel_vals = VT;
[p, fiteval] = polyfit(VTsel_vals, Vsel_vals, 1);
rho = corr(VTsel_vals, Vsel_vals, "Type","Spearman");
xfit = [min(VTsel_vals), max(VTsel_vals)];
yfit = polyval(p, xfit);
plot(xfit, yfit, mkr{2,v-1}, 'Color', clr{s}.^g, 'LineWidth', 1.5);
txt1 = ['  \it y\rm = ',num2str(p(1),'%.1f'),'\it x\rm + ',num2str(p(2),'%.1f'),'  '];
%txt2 = ['  R^2 = ',num2str(fiteval.rsquared,'%+.2f'),'  '];
txt2 = ['  \rho = ',num2str(rho,'%+.2f'),'  '];
if alignR
    text(xfit(end), yfit(end), txt1, ...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
        'Color',clr{s}.^g);
    text(xfit(end), yfit(end), txt2, ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        'Color',clr{s}.^g);
else
    text(xfit(1), yfit(1), txt1, ...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom', ...
        'Color',clr{s}.^g);
    text(xfit(1), yfit(1), txt2, ...
        'HorizontalAlignment','right', 'VerticalAlignment','top', ...
        'Color',clr{s}.^g);
end
alignR = ~alignR;

            end
        end
    end

grid on;
title("Subject: "+subj);
xlabel('\it x\rm = XTRA VT'); ylabel(['\it y\rm = ',bandname]);
legend(lgd, 'Location','eastoutside')
xl = xlim; 
xl = xl + [-1,1]*diff(xl)*0.25;
xlim(xl);
end