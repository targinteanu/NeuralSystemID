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

%% batch run and plot 

tblsAll = cell(size(fn));

for subjtype = 1:length(fn)
    fn_ = fn{subjtype};
    tblsThisGroup = [repmat({table}, 2, 2); ...
                     cell(1, 2)];
        % rows = {Avg; Std; Num} 
        % cols = {source, sink}

    for subj = 1:length(fn_) 
        fn_subj = fn_{subj};
        load([fp,filesep,fn_subj]);
        
        subj_name = fn_subj; 
        subj_name_end = strfind(subj_name, ' --- ');
        if ~isempty(subj_name_end)
            subj_name_end = subj_name_end(1); 
            subj_name = subj_name(1:subj_name_end);
        end
        clear subj_name_end
        disp(['Processing ',subj_name])

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
    tblsAll{subjtype} = tblsThisGroup;
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