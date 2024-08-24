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

%% all channels and subjects
trainEval = {}; testEval = {};

for fn_ = fn
    fn_ = fn_{:};
    sn_ = {};
    for fn_subj = fn_
        fn_subj = fn_subj{:};
        load([fp,filesep,fn_subj], 'EEG_table');
        eeg = EEG_table.BaselineOpen('before experiment'); 

% epoch 
[EEGlist{1}, EpoList{1}] = epochBaseline(EEG_table,'BaselineOpen','both',1,1); 
[EEGlist{2}, EpoList{2}] = epochBaseline(EEG_table,'BaselineClosed','both',1,1); 
[EEGlist{3}, EpoList{3}] = epochStim(EEG_table,'TempStim',1);
[EEGlist{4}, EpoList{4}] = epochStim(EEG_table,'PinPrick',1);

% fit 
trnE = cell(size(EEGlist)); tstE = cell(size(EEGlist));
for i = 1:length(EEGlist)
[A,trnE{i},tstE{i}] = ...
    plotModelFit(EEGlist{i}, EpoList{i}, ...
    @(tsTbl, trTbl) fitLTIauton(tsTbl, trTbl), ...
    '', [13,30]);
end
trainEval = [trainEval; trnE]; testEval = [testEval; tstE];

    end
end

%% analysis 
tstEvalMean = cellfun(@(x)mean([x.pRMSE]), testEval); 
trnEvalMean = cellfun(@(x)mean([x.pRMSE]), trainEval); 
evalMean = {tstEvalMean; trnEvalMean};
figure;
for m = 1:length(evalMean)
    TE = evalMean{m};
    subplot(2,1,m);
    for i = 1:size(TE,2)
        histogram(TE(:,i)); hold on;
    end
    legend('BaselineOpen', 'BaselineClosed', 'TempStim', 'PinPrick')
    xlabel('mean pRMSE')
end
subplot(2,1,1); title('test');
subplot(2,1,2); title('train');