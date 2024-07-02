%% get EEG list to load

basedir = cd; 
cd('/Users/torenarginteanu/Desktop/Data_Chronic Pain'); 
fp = uigetdir('.mat', 'Select Pre- or Post-processed EEG folder');
cd(basedir); 

files = dir(fp);
files = files(~[files.isdir]);
fn = {files.name};

fnH = listdlg("PromptString",'Select Controls', ...
    'ListString',fn, 'SelectionMode','multiple');
fnP = listdlg("PromptString",'Select Patients', ...
    'ListString',fn, 'SelectionMode','multiple');
fnH = fn(fnH); fnP = fn(fnP); 

fn = {fnH, fnP}; 
groupname = {'Control', 'Patient'};

%% get list of all channels and subjects
chlbl_all = {}; sn = {};
chloc_all = [];
for fn_ = fn
    fn_ = fn_{:};
    sn_ = {};
    for fn_subj = fn_
        fn_subj = fn_subj{:};
        load([fp,filesep,fn_subj], 'EEG_table');
        eeg = EEG_table.BaselineOpen('before experiment'); 
        eeg = eeg{1}(1); 
        chlbl_all = [chlbl_all, upper({eeg.chanlocs.labels})];
        chloc_all = [chloc_all, eeg.chanlocs];

        subj_name = fn_subj;
        subj_name_end = strfind(subj_name, ' --- ');
        if ~isempty(subj_name_end)
            subj_name_end = subj_name_end(1)-1; 
            subj_name = subj_name(1:subj_name_end);
        end
        sn_ = [sn_, subj_name];
        clear subj_name_end subj_name
    end
    sn = [sn, {sn_}];
end
[chlbl_all, ui] = unique(chlbl_all); chloc_all = chloc_all(ui);
clear fn_ fn_subj EEG_table eeg sn_ ui

%% batch computations 

tblsAll = cell(size(fn));

for subjgroup = 1:length(fn)
    fn_ = fn{subjgroup};
    sn_ = sn{subjgroup};

    blanktable = table('Size', [length(fn_), length(chlbl_all)], ...
        'VariableTypes', repmat("single",size(chlbl_all)), ...
        'VariableNames', chlbl_all, 'RowNames', sn_);
    blanktable{:,:} = nan;
    tblsThisGroup = [repmat({blanktable}, 2, 2); ...
                     cell(1, 2)];
        % rows = {Avg; Std; Num} 
        % cols = {source, sink}

    for subj = 1:length(fn_) 
        fn_subj = fn_{subj};
        subj_name = sn_{subj};
        disp(['Processing ',subj_name])
        load([fp,filesep,fn_subj], 'EEG_table'); 

%% epoch 
%[curEEGlist, EpocList] = epochStim(EEG_table,'TempStim',1);
[curEEGlist, EpocList] = epochBaseline(EEG_table,'BaselineClosed','both',1,1);

%% fit and source-sink 

A = plotModelFit(curEEGlist, EpocList, ...
    @(tsTbl, trTbl) fitLTIauton(tsTbl, trTbl), ...
    '');
chanlocs = curEEGlist(1).chanlocs;
[srcness,snkness] = SourceSink(A, chanlocs, false);
chlbl = upper({chanlocs.labels});

srcAvg = mean(srcness,3); 
srcStd = std(srcness,[],3); 
srcNum = size(srcness,3); 

snkAvg = mean(snkness,3)'; 
snkStd = std(snkness,[],3)'; 
snkNum = size(snkness,3)'; 

tblsThisGroup{1,1} = addToTable(tblsThisGroup{1,1}, ...
    {subj_name}, chlbl, srcAvg);
tblsThisGroup{2,1} = addToTable(tblsThisGroup{2,1}, ...
    {subj_name}, chlbl, srcStd);
tblsThisGroup{3,1} = [tblsThisGroup{3,1}; srcNum];

tblsThisGroup{1,2} = addToTable(tblsThisGroup{1,2}, ...
    {subj_name}, chlbl, snkAvg);
tblsThisGroup{2,2} = addToTable(tblsThisGroup{2,2}, ...
    {subj_name}, chlbl, snkStd);
tblsThisGroup{3,2} = [tblsThisGroup{3,2}; snkNum];

    end
    tblsAll{subjgroup} = tblsThisGroup;
end

%% cumulative results 

for subjgroup = 1:length(fn)
    tblsThisGroup = tblsAll{subjgroup};
    for m = 1:size(tblsThisGroup,2) % {source, sink}
        % 1) Avg -- mean 
        tbl = tblsThisGroup{1,m}; 
        tbl = addToTable(tbl, {'Cumulative'}, tbl.Properties.VariableNames, ...
            mean(tbl{:,:})); % omitnan? 
        tblsThisGroup{1,m} = tbl;
        % 2) Std -- rms
        tbl = tblsThisGroup{2,m}; 
        tbl = addToTable(tbl, {'Cumulative'}, tbl.Properties.VariableNames, ...
            rms(tbl{:,:})); % omitnan? 
        tblsThisGroup{2,m} = tbl;
        % 3) Num -- sum 
        tblsThisGroup{3,m} = [tblsThisGroup{3,m}; sum(tblsThisGroup{3,m})];
    end
    tblsAll{subjgroup} = tblsThisGroup;
end

%% plot results 

srcfig = figure('Units','normalized', 'Position',[.05,.05,.9,.4]);
sgtitle('Source-ness');
snkfig = figure('Units','normalized', 'Position',[.05,.5,.9,.4]);
sgtitle('Sink-ness');

% 4 rows: mean-H, std or SE, mean-P, std or SE
% last column: cumulative H/P
fn_len = cellfun(@(lst) length(lst), fn);
W = max(fn_len) + 1; H = 4;

for subjgroup = 1:length(fn)
    tblsThisGroup = tblsAll{subjgroup};
    for m = 1:size(tblsThisGroup,2) % {source, sink}

        if m==1
            figure(srcfig); 
        else
            figure(snkfig);
        end

        tblAvg = tblsThisGroup{1,m};
        tblStd = tblsThisGroup{2,m};
        N = tblsThisGroup{3,m};

        [tblAvg,keptchan] = pruneNan(tblAvg); 
        tblStd = pruneNan(tblStd); 
        chloc = chloc_all(keptchan);

        sn_ = tblAvg.Properties.RowNames;
        for subj = 1:height(tblAvg)
            subj_name = sn_{subj};
            subj_name = subj_name(1:4); % or shorten to first space?
            r1 = W*2*(subjgroup-1); r2 = r1+W; 
            subplot(H,W,r1+subj); 
            topoplot(tblAvg{subj,:}, chloc, 'maplimits', 'maxmin'); 
            title(subj_name); colorbar; 
            subplot(H,W,r2+subj); 
            topoplot(tblStd{subj,:}, chloc, 'maplimits', 'maxmin');
            title(subj_name); colorbar; 
        end
    end
end

%% helpers 

function tbl = addToTable(tbl, rowNames, colNames, data)
for r = 1:size(data,1)
    rowName = rowNames{r};
    rowName = ['''',rowName,''''];
    for c = 1:size(data,2)
        colName = colNames{c};
        data_ = data(r,c);
        try
        eval(['tbl.',colName,'(',rowName,') = data_;']);
        catch ME
            ME.message
        end
    end
end
end

function [tbl, keptcols] = pruneNan(tbl)
data = tbl{:,:}; 
keptcols = sum(isnan(data)); 
keptcols = ~(keptcols); 
tbl = tbl(:,keptcols); 
end