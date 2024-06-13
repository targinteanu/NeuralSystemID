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
for fn_ = fn
    fn_ = fn_{:};
    sn_ = {};
    for fn_subj = fn_
        fn_subj = fn_subj{:};
        load([fp,filesep,fn_subj], 'EEG_table');
        eeg = EEG_table.BaselineOpen('before experiment'); 
        eeg = eeg{1}(1); 
        chlbl_all = [chlbl_all, upper({eeg.chanlocs.labels})];

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
clear fn_ fn_subj EEG_table eeg sn_
chlbl_all = unique(chlbl_all);

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
        load([fp,filesep,fn_subj]); 

%% epoch 

curEEGlist = [EEG_table.PinPrick('before experiment'), EEG_table.PinPrick('after experiment')];
curEEGlist = [curEEGlist{:}]; 

    % Determine the epoch duration and overlap: 
    epochT = 1; % s
    epochT = .5*epochT;

EpocList = cell(size(curEEGlist));
for lstIdx = 1:length(curEEGlist)
    eeg = curEEGlist(lstIdx);
    if ~isempty(eeg)
        curEpochs3D = pop_epoch(eeg, {'11'}, [-epochT, epochT]);
        t = 1:curEpochs3D.trials; 
        curEpochs = repmat(curEpochs3D, size(t));
        for idx = 1:length(t)
            curEpoch = curEpochs(idx); 
            curEpoch.trials = 1;
            curEpoch.data = curEpoch.data(:,:,idx);
            curEpochs(idx) = curEpoch;
        end
        EpocList{lstIdx} = curEpochs;
    end
end

clear curEpoch curEpochs eeg idx lstIdx t 

%% fit and source-sink 

A = plotModelFit(curEEGlist, EpocList, ...
    @(tsTbl, trTbl) fitLTIauton(tsTbl, trTbl), ...
    '');
chanlocs = curEEGlist(1).chanlocs;
[srcness,snkness,figs] = SourceSink(A, chanlocs, false);
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
% to do: open all tables in tblsAll 
% make a cumulative row 
% for Avg: use mean 
% for Std: use rms 
% for Num: use sum 
% omitnan? 
% prune nan columns of all tables? 

%% plot results 
% to do: open all tables in tblsAll 
% 4 rows: mean-H, std or SE, mean-P, std or SE
% last column: cumulative H/P
% title with subj name or "cumulative"
% shorten subj name further? (to first space)

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