%% setup

% find files and determine data fields 
fp = uigetdir;
filelist = dir(fullfile(fp, '*.mat'));
file1 = filelist(~[filelist.isdir]);
file1 = file1(arrayfun(@(f) f.bytes > 0, file1));
file1 = file1(arrayfun(@(f) f.name(1)~='.', file1));

channelNamesAll = []; Fs = [];
for fi = 1:length(file1)

curfile = file1(fi); 
curfiledata = load(fullfile(curfile.folder, curfile.name));
channelNames = getChannelNames(curfiledata);

% channel selection 
% For now, only select ANALOG_IN and LFP channels; ignore spike, RAW, and
% SEG. In the future, spike sorting should be performed later in this
% script instead of throwing out spikes. Not sure what to do with RAW/SEG
chincl = contains(channelNames, 'LFP') | ...
         contains(channelNames, 'ANALOG_IN');
channelNames = channelNames(chincl); 

% get sample rates
fs = cellfun(@(chN) getsamplerate(chN, curfiledata), channelNames);

% store, clear, move on
Fs = [Fs; fs];
channelNamesAll = [channelNamesAll; channelNames];
clear curfile curfiledata
clear channelNames fs

end
channelNames = channelNamesAll; clear channelNamesAll

[channelNames, iFs] = unique(channelNames);
Fs = Fs(iFs); clear iFs
channelNames = channelNames'; Fs = Fs';

% group channels by sampling rate 
FsGrouped = unique(Fs);
channelNamesGrouped = cell(size(FsGrouped));
for grp = 1:length(FsGrouped)
    grpInds = Fs == FsGrouped(grp);
    channelNamesGrouped{grp} = channelNames(grpInds);
end
clear grp grpInds

fileDataFields = {...
    'SF_DRIVE_CONF', ...
    ... 'SF_DRIVE_DOWN', ...
    'SF_DRIVE_SET_POS', ...
    ... 'SF_DRIVE_STOP', ...
    ... 'SF_DRIVE_UP', ...
    'SF_LEVEL' ...
    }; 
channelDataFields = {'BitResolution', 'Gain'};

Tbls0 = cell(size(FsGrouped));
Tbls = Tbls0;

% set file saving location
svloc = [fp,filesep,'Saved To Table',filesep,'Table Data ',...
    datestr(datetime, 'yyyy-mm-dd HH.MM.SS'),' ',...
    getFileVersion(mfilename("fullpath"))];
svN = 1;
sizethresh = 2e10; % size (bytes) at which to save and clear
pause(1)
mkdir(svloc); 

%% run1 - populate "horizontal" 

for f = filelist'
    %{
    clearvars -except ...
        Tbls Tbls0 channelNames channelDataFields fileDataFields ...
        f filelist svloc svN sizethresh ...
        Fs FsGrouped channelNamesGrouped
    %}
    if (~f.isdir) && (f.bytes > 0)
        fnfull = fullfile(f.folder, f.name); 
        [fp,fn,fe] = fileparts(fnfull);
        fni = find(fn == '.');
        fn1 = fn(1:(fni-1)); fn2 = fn((fni+1):end);
        fn1 = string(fn1); fn2 = string(fn2);
        if strcmpi(fe, '.mat')
            curfiledata = load(fnfull);

            %{
            % anticipate if the memory will become full; 
            % if so, save and clear 
            sz = whos('Tbls'); sz = sz.bytes;
            if (sz + f.bytes > sizethresh) && (f.bytes < sizethresh)
                save([svloc,filesep,'SavedTables',num2str(svN),'.mat'], "Tbls","channelNames");
                svN = svN+1;
                clear Tbls
                Tbls = Tbls0;
            end
            %}

            % look at channels 
            chNames = getChannelNames(curfiledata);
            chincl = contains(chNames, 'LFP') | ...
                     contains(chNames, 'ANALOG_IN');
            chNames = chNames(chincl);
            if ~strcmpwrapper(channelNames, chNames)
                warning(['On ',fn,': channel names do not match'])
            end
            FileData = varnames2struct(fileDataFields, curfiledata, '');
            for FSGRP = 1:length(channelNamesGrouped)
                TT = []; T1 = inf; T2 = -inf;
                chNamesGrp = channelNamesGrouped{FSGRP};
                for CHNAME = chNamesGrp
                chName = CHNAME{:};
                try
                t1 = curfiledata.([chName,'_TimeBegin']); % s
                t2 = curfiledata.([chName,'_TimeEnd']); % s
                fs = curfiledata.([chName,'_KHz'])*1000; % Hz
                T1 = min(T1, t1); T2 = max(T2, t2);
                t = t1:(1/fs):t2; 
                if isfield(curfiledata, [chName,'_BitResolution'])
                    DataRes = curfiledata.([chName,'_BitResolution']);
                else
                    DataRes = 1;
                end
                if isfield(curfiledata, [chName,'_Gain'])
                    DataGain = curfiledata.([chName,'_Gain']);
                else
                    DataGain = 1;
                end
                Data = curfiledata.(chName); % int
                Data = double(Data)*DataRes/DataGain;
                if length(t) ~= length(Data)
                    err = abs(length(Data) - length(t))/length(Data);
                    warning([chName,' mismatch in time and data length by ', num2str(100*err),'%'])
                    % treat times as more accurate than sample rate 
                    t = linspace(t1, t2, length(Data));
                end

                t = seconds(t);
                T = array2timetable(Data', 'RowTimes', t'); T.Properties.VariableNames{1} = chName;
                if isempty(TT)
                    TT = T;
                else
                    TT = synchronize(TT, T);
                end
                
                catch ME
                    warning(['On ',fn,' - ',chName,': ',ME.message])
                end
                clear chName t1 t2 t Data err FileName T
                end

                % mark table with file details 
                TT.Properties.Events = eventtable(seconds([T1;T2]), ...
                    "EventLabels", string(fn)+[" Start";" End"]);

                Tbls{FSGRP} = tblvertcat(Tbls{FSGRP}, TT);
                % if the memory is getting full, save and clear 
                sz = whos('Tbls'); sz = sz.bytes;
                if sz > sizethresh
                    save([svloc,filesep,'SavedTables',num2str(svN),'.mat'], "Tbls");
                    svN = svN+1;
                    clear Tbls
                    Tbls = Tbls0;
                end
            end
            clear curfiledata FileData
        end
    end
end

% save if haven't yet
if sum(cellfun(@numel, Tbls))
    save([svloc,filesep,'SavedTables',num2str(svN),'.mat'], "Tbls");
    svN = svN+1;
end
clear Tbls
Tbls = Tbls0;

%% run2 - consolidate saved files "vertically" 
% if svN > 2
% there are multiple "Horizontal" files to consolidate

clearvars -except channelNames svloc svN sizethresh ...
    FsGrouped channelNamesGrouped
channelNames0 = channelNames; clear channelNames
filelist = dir(svloc);
filelist(~[filelist.isdir]);
filelist = filelist(arrayfun(@(f) f.bytes > 0, filelist));
filelist = filelist(contains({filelist.name}, 'SavedTables'));
svN = 1;

for FSGRP = 1:length(channelNamesGrouped)
    fs = FsGrouped(FSGRP);
    Tbl = [];
    for f = filelist'
        fnfull = fullfile(f.folder, f.name); 
        [fp,fn,fe] = fileparts(fnfull);
        if strcmpi(fe, '.mat')
            load(fnfull);
            T = Tbls{FSGRP}; clear Tbls;
            if isempty(Tbl)
                Tbl = T;
            else
                Tbl = tblvertcat(Tbl, T);
            end

            % if the memory is getting full, save and clear
            sz = whos('Tbl'); sz = sz.bytes;
            if sz > sizethresh
                Tbl = sortrows(Tbl, 'Time');
                save([svloc,filesep,'SavedTable',num2str(fs),'Hz',num2str(svN),'.mat'], "Tbl");
                svN = svN+1;
                clear Tbl
                Tbl = [];
            end
        end
    end

    % save if haven't yet 
    if numel(Tbl)
        Tbl = sortrows(Tbl, 'Time');
        save([svloc,filesep,'SavedTable',num2str(fs),'Hz',num2str(svN),'.mat'], "Tbl");
        svN = svN+1;
    end
end

% end

%% helper 

function fs = getsamplerate(chName, filedata)
fs = 0;
    try 
        fs = filedata.([chName,'_KHz'])*1000; % Hz
    catch ME
        warning(['On ',chName,': ',ME.message])
    end
end

function channelNames = getChannelNames(curfiledata)
if isfield(curfiledata, 'Channel_ID_Name_Map')
    channelNames = {curfiledata.Channel_ID_Name_Map.Name};
else
    allFields = fieldnames(curfiledata);
    channelFields = cellfun(@(str) str(1)=='C', allFields);
    channelNames = allFields(channelFields);
    channelInfo = ... trim away fields that are info about channels
        contains(lower(channelNames), 'hz') | ...
        contains(lower(channelNames), 'resolution') | ...
        (contains(lower(channelNames), 'level') & ~contains(lower(channelNames), 'level_seg')) | ...
        contains(lower(channelNames), 'gain') | ...
        contains(lower(channelNames), 'time');
    channelNames = channelNames(~channelInfo);
    % TO DO: 'LEVEL_SEG' should be removed from the end of some channel names!
end
end

function S = varnames2struct(varnames, filedata, header)
if nargin < 3
    header = '';
end
vars = cellfun(@(vn) getfieldwrapper(filedata, [header,vn]), varnames, 'UniformOutput', false);
S = cell2struct(vars, varnames, 2);
end

function fld = getfieldwrapper(S, fldname)
if isfield(S, fldname)
    fld = S.(fldname);
else
    warning(['Unrecognized field name ',fldname]);
    fld = [];
end
end

function yn = strcmpwrapper(strs1, strs2)
yn = length(strs1) == length(strs2); 
if yn 
    yn = prod(strcmp(strs1, strs2));
end
end

function T = tblvertcat(T1, T2)
try
    T = [T1; T2];
catch ME
    if contains(ME.identifier, 'SizeMismatch')
        T1small = width(T1) < width(T2);
        if T1small
            Tsmall = T1; Tlarge = T2;
        else
            Tsmall = T2; Tlarge = T1;
        end
        smallvars = Tsmall.Properties.VariableNames;
        for VAR = Tlarge.Properties.VariableNames
            var = VAR{:};
            if ~sum(strcmp(var, smallvars))
                Tsmall.(var) = nan(height(Tsmall),1);
            end
        end
        if T1small
            T = tblvertcat(Tsmall, Tlarge);
        else
            T = tblvertcat(Tlarge, Tsmall);
        end
    else
        rethrow(ME)
    end
end
end