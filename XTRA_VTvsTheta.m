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
for s = 1:length(subjnames)
    subj = subjnames(s);
    VT = data2{s,1};
    if ~isempty(VT)
        BL1 = data2{s,2};
        NB = data2{s,3};
        BL2 = data2{s,4};
        if ~isempty(BL1)
            plot(VT, BL1, '.r'); hold on; 
        end
        if ~isempty(NB)
            plot(VT, NB, '.b'); hold on;
        end 
        if ~isempty(BL2)
            plot(VT, BL2, '.m'); hold on;
        end
    end
end
grid on;
xlabel('XTRA VT'); ylabel('Theta Power (dB)');