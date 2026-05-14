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
            elseif strcmp(vname, 'NB')
                vi = 3;
            elseif strcmp(vname, 'BL2')
                vi = 4;
            else
                warning(['Unrecognized variable name ',vname]);
            end
            data2{s, vi} = [data{3:end, v}]';
        end
    end
end

figure; 
mkr = {'s', 'x', 'o'};
clr = {[0.0660    0.4430    0.7450], ... blue 
       [0.8660    0.3290         0], ... red 
       [0.2310    0.6660    0.1960], ... green
       [0.5210    0.0860    0.8190], ... purple
       [0.4645    0.3470    0.0625], ... dark yellow
       ...[0.9290    0.6940    0.1250], ... yellow
       [0.1840    0.7450    0.9370], ... teal
       [0.8190    0.0150    0.5450]  ... dark red
       };
for s = 1:length(subjnames)
    subj = subjnames(s);
    VT = data2{s,1};
    if ~isempty(VT)
        for v = 2:4
            V = data2{s,v};
            if ~isempty(V)
                plot(VT, V, '.', ...
                    'Marker', mkr{v-1}, 'Color', clr{s}); 
                hold on;
            end
        end
    end
end
grid on;
xlabel('XTRA VT'); ylabel('Theta Power (dB)');