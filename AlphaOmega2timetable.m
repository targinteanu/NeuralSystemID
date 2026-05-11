%% setup

% find files and determine data fields 
fp = uigetdir;
filelist = dir(fullfile(fp, '*.mat'));
file1 = filelist(~[filelist.isdir]);
file1 = file1(arrayfun(@(f) f.bytes > 0, file1));
file1 = file1(arrayfun(@(f) f.name(1)~='.', file1));
%%
channelNames = []; Fs = []; channelIDs = [];
stimNames = []; stimFs = []; stimIDs = [];
mindepth = inf; maxdepth = -inf;
mintime = inf; maxtime = -inf;
% stim IDs appear to be independent of stim channel IDs; unclear if it is
% necessary to track these or stimFs
for fi = 1:length(file1)

curfile = file1(fi); 
curfiledata = load(fullfile(curfile.folder, curfile.name));
[chNames, chIDs] = getChannelNames(curfiledata);
[stNames, stIDs] = getStimNames(curfiledata);

% get details of this rec 
[data, flag] = sscanf(curfile.name, '%ct%fd%ff%f %s');
if flag > 4
    SIDE = string(char(data(1)));
    N = data(2);
    DEPTH = data(3);
    FILE = data(4);
    mindepth = min(mindepth, DEPTH);
    maxdepth = max(maxdepth, DEPTH);
end

% get sample rates
fs = cellfun(@(chN) getsamplerate(chN, curfiledata), chNames);
stimfs = cellfun(@(chN) getsamplerate(chN, curfiledata), stNames);

% get times 
for chName = chNames
    chname = chName{:};
    if isfield(curfiledata, [chname,'_TimeBegin'])
        mintime = min(mintime, curfiledata.([chname,'_TimeBegin']));
    end
    if isfield(curfiledata, [chname,'_TimeEnd'])
        maxtime = max(maxtime, curfiledata.([chname,'_TimeEnd']));
    end
end

% store, clear, move on
Fs = [Fs, fs];
channelNames = [channelNames, chNames];
channelIDs = [channelIDs, chIDs];
stimFs = [stimFs, stimfs];
stimNames = [stimNames, stNames];
stimIDs = [stimIDs, stIDs];
clear curfile curfiledata
clear chNames chIDs fs stNames stIDs stimfs
end

[channelNames, iCh] = unique(channelNames);
Fs = Fs(iCh); channelIDs = channelIDs(iCh); clear iCh

[stimNames, iSt] = unique(stimNames);
stimFs = stimFs(iSt); stimIDs = stimIDs(iSt); clear iSt

% channel selection manually
channelNamesWithFs = arrayfun(@(i) [channelNames{i},': ',num2str(Fs(i)),'Hz'], ...
    1:length(Fs), 'UniformOutput',false);
[chincl, chselmade] = listdlg("PromptString",'Select Data Channel(s)', ...
    "SelectionMode","multiple", ...
    "ListString",channelNamesWithFs, "ListSize",[300 500]); 
if ~chselmade
    error('Selection must be made. Use Select All to choose all channels.')
end
channelNames = channelNames(chincl); Fs = Fs(chincl);
channelIDs = channelIDs(chincl);
channelNamesWithFs = channelNamesWithFs(chincl); % used anywhere else?

% depth selection manually 
minmaxdepth = inputdlg({'Minimum Depth:', 'Maximum Depth:'}, ...
    'Specify Depth Range', 1, {num2str(mindepth), num2str(maxdepth)});
if ~isempty(minmaxdepth)
    minDepth = str2double(minmaxdepth{1});
    maxDepth = str2double(minmaxdepth{2});
end
% time selection manually 
minmaxtime = inputdlg({'Start Time (s):', 'End Time (s):'}, ...
    'Specify Time Range', 1, {num2str(mintime), num2str(maxtime)});
if ~isempty(minmaxtime)
    mintime = str2double(minmaxtime{1});
    maxtime = str2double(minmaxtime{2});
end

% group channels by sampling rate 
FsGrouped = unique(Fs);
channelNamesGrouped = cell(size(FsGrouped, 2)); % [name, ID] 
for grp = 1:length(FsGrouped)
    grpInds = Fs == FsGrouped(grp);
    channelNamesGrouped{grp,1} = channelNames(grpInds);
    channelNamesGrouped{grp,2} = channelIDs(grpInds);
end
clear grp grpInds minmaxdepth minmaxtime
mintime = seconds(mintime); maxtime = seconds(maxtime);

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
            FileData = varnames2struct(fileDataFields, curfiledata, '');

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

            % get details of this rec
            [data, flag] = sscanf(fn, '%ct%fd%ff%f %s');
            if flag > 4
                SIDE = string(char(data(1)));
                N = data(2);
                DEPTH = data(3);
                FILE = data(4);
            end

            if (DEPTH >= mindepth) && (DEPTH <= maxdepth)

            % process stim markers into event table 
            ET = [];
            for STNAME = stimNames
                try
                stName = STNAME{:};
                stimTimesIDs = curfiledata.(stName);
                stimfs = curfiledata.([stName,'_KHz'])*1000; % Hz
                stimTimes = stimTimesIDs(1,:) / stimfs; % s
                chID = stimTimesIDs(2,:);
                chNames = cell(size(chID));
                for iSt = 1:length(chID)
                    iCh = channelIDs == chID(iSt);
                    if sum(iCh) == 1
                        chName = channelNames{find(iCh)};
                    else
                        chName = num2str(chID(iSt));
                    end
                    chNames{iSt} = [stName,' on channel ',chName];
                end
                stimTimes = stimTimes'; chNames = chNames';
                E = eventtable(seconds(stimTimes), "EventLabels",string(chNames));
                if isempty(ET)
                    ET = E;
                else
                    ET = [ET; E];
                end
                clear E stName stimTimesIDs sitmTimes stimfs chID chNames 
                clear iSt iCh chID chName
                catch ME
                    warning(ME.message)
                end
            end

            %{
            % look at channels 
            if ~strcmpwrapper(channelNames, chNames)
                warning(['On ',fn,': channel names do not match'])
            end
            %}
            
            % process a group of channels with same Fs into TT
            for FSGRP = 1:size(channelNamesGrouped,1)
                TT = []; T1 = inf; T2 = -inf;
                chNamesGrp = channelNamesGrouped{FSGRP,1};
                % process a single channel in the group into T
                for CHNAME = chNamesGrp
                chName = CHNAME{:};
                %if sum(strcmp(channelNames, chName))
                try

                % interpret data
                % sometimes there is no TimeBegin/End - why?? 
                if isfield(curfiledata, chName) && (...
                        ~isfield(curfiledata,[chName,'_TimeBegin']) || ...
                        ~isfield(curfiledata,[chName,'_TimeEnd']) || ...
                        ~isfield(curfiledata,[chName,'_KHz']))
                    keyboard
                end
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

                % construct T and append it to TT 
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
                %end
                clear chName t1 t2 t Data err FileName T
                end

                % mark table TT with file details 
                TT.Properties.Events = eventtable(seconds([T1;T2]), ...
                    "EventLabels", string(fn)+[" Start";" End"]);

                % include stims as events 
                if ~isempty(ET)
                    TT.Properties.Events = [TT.Properties.Events; ET];
                end

                % limit time window 
                timesel = (TT.Time >= mintime) & (TT.Time <= maxtime);
                TT = TT(timesel,:);

                if ~isempty(TT)
                    Tbls{FSGRP} = tblvertcat(Tbls{FSGRP}, TT);
                end
                % if the memory is getting full, save and clear 
                sz = whos('Tbls'); sz = sz.bytes;
                if sz > sizethresh
                    disp(['Saving ',svloc,filesep,'SavedTables',num2str(svN),'.mat'])
                    save([svloc,filesep,'SavedTables',num2str(svN),'.mat'], "Tbls");
                    svN = svN+1;
                    clear Tbls
                    Tbls = Tbls0;
                end

            end
            end
            clear curfiledata FileData
        end
    end
end

% save if haven't yet
if sum(cellfun(@numel, Tbls))
    disp(['Saving ',svloc,filesep,'SavedTables',num2str(svN),'.mat'])
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

for FSGRP = 1:size(channelNamesGrouped,1)
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
                disp(['Saving ',svloc,filesep,'SavedTable',num2str(fs),'Hz',num2str(svN),'.mat'])
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
        disp(['Saving ',svloc,filesep,'SavedTable',num2str(fs),'Hz',num2str(svN),'.mat'])
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

function [channelNames, channelIDs] = getChannelNames(curfiledata)
if isfield(curfiledata, 'Channel_ID_Name_Map')
    channelNames = {curfiledata.Channel_ID_Name_Map.Name};
    channelIDs = [curfiledata.Channel_ID_Name_Map.ID];
else
    allFields = fieldnames(curfiledata);
    channelFields = cellfun(@(str) str(1)=='C', allFields);
    channelNames = allFields(channelFields);
    channelInfo = ... trim away fields that are info about channels
        contains(lower(channelNames), 'hz') | ...
        contains(lower(channelNames), 'stimmarker') | ...
        contains(lower(channelNames), 'resolution') | ...
        (contains(lower(channelNames), 'level') & ~contains(lower(channelNames), 'level_seg')) | ...
        contains(lower(channelNames), 'gain') | ...
        contains(lower(channelNames), 'time');
    channelNames = channelNames(~channelInfo);
    channelNames = channelNames'; % horizontal
    channelIDs = nan(size(channelNames));
    % TO DO: 'LEVEL_SEG' should be removed from the end of some channel names!
end
end

function [stimNames, stimIDs] = getStimNames(curfiledata)
if isfield(curfiledata, 'Ports_ID_Name_Map')
    allNames = {curfiledata.Ports_ID_Name_Map.Name};
    allIDs = [curfiledata.Ports_ID_Name_Map.ID];
    stimFields = contains(allNames, 'CStimMarker');
    stimNames = allNames(stimFields);
    stimIDs = allIDs(stimFields);
else
    allFields = fieldnames(curfiledata);
    stimFields = contains(allFields, 'CStimMarker');
    stimNames = allFields(stimFields);
    channelInfo = ... trim away fields that are info about channels
        contains(lower(stimNames), 'hz') | ...
        contains(lower(stimNames), 'resolution') | ...
        (contains(lower(stimNames), 'level') & ~contains(lower(stimNames), 'level_seg')) | ...
        contains(lower(stimNames), 'gain') | ...
        contains(lower(stimNames), 'time');
    stimNames = stimNames(~channelInfo);
    stimNames = stimNames'; % horizontal
    stimIDs = nan(size(stimNames)); % to be filled in later
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
    yn = prod(strcmp(sort(strs1), sort(strs2)));
end
end

function T = tblvertcat(T1, T2)
try
    T = [T1; T2];
catch ME
    if contains(ME.identifier, 'SizeMismatch') || contains(ME.identifier, 'UnequalVarNames')
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