%% setup

% find files and determine data fields 
%fp = uigetdir;
filelist = dir(fp);
file1 = filelist(~[filelist.isdir]);
file1 = file1(arrayfun(@(f) f.bytes > 0, file1));
file1 = file1(arrayfun(@(f) f.name(1)~='.', file1));
file1 = file1(1); load(fullfile(file1.folder, file1.name));
channelNames = {Channel_ID_Name_Map.Name};
Tbls0 = cell(size(channelNames));
Tbls = Tbls0;
fileDataFields = {...
    'SF_DRIVE_CONF', ...
    ... 'SF_DRIVE_DOWN', ...
    'SF_DRIVE_SET_POS', ...
    ... 'SF_DRIVE_STOP', ...
    ... 'SF_DRIVE_UP', ...
    'SF_LEVEL' ...
    }; 
channelDataFields = {'BitResolution', 'Gain'};

% channel selection 
% For now, only select ANALOG_IN and LFP channels; ignore spike, RAW, and
% SEG. In the future, spike sorting should be performed later in this
% script instead of throwing out spikes. Not sure what to do with RAW/SEG
chincl = contains(channelNames, 'LFP') | ...
         contains(channelNames, 'ANALOG_IN');
channelNames = channelNames(chincl); 

% group channels by sampling rate 
Fs = cellfun(@getsamplerate, channelNames);
FsGrouped = unique(Fs);
channelNamesGrouped = cell(size(FsGrouped));
for grp = 1:length(FsGrouped)
    grpInds = Fs == FsGrouped(grp);
    channelNamesGrouped{grp} = channelNames(grpInds);
end
clear grp grpInds

% set file saving location
svloc = [fp,filesep,'Saved To Table',filesep,'Table Data ',...
    datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];
svN = 1;
sizethresh = 1e9; % size (bytes) at which to save and clear
pause(1)
mkdir(svloc); 

%% run1 - populate "horizontal" 

for f = filelist'
    clearvars -except ...
        Tbls Tbls0 channelNames channelDataFields fileDataFields ...
        f filelist svloc svN sizethresh ...
        Fs FsGrouped channelNamesGrouped
    if (~f.isdir) && (f.bytes > 0)
        fnfull = fullfile(f.folder, f.name); 
        [fp,fn,fe] = fileparts(fnfull);
        fni = find(fn == '.');
        fn1 = fn(1:(fni-1)); fn2 = fn((fni+1):end);
        fn1 = string(fn1); fn2 = string(fn2);
        if strcmpi(fe, '.mat')
            load(fnfull)

            %{
            % anticipate if the memory will become full; 
            % if so, save and clear 
            sz = whos('Tbls'); sz = sz.bytes;
            if (sz + f.bytes > sizethresh) && (f.bytes < sizethresh)
                save([svloc,filesep,'SaveFileH',num2str(svN),'.mat'], "Tbls","channelNames");
                svN = svN+1;
                clear Tbls
                Tbls = Tbls0;
            end
            %}

            % look at channels 
            if ~strcmpwrapper(channelNames, {Channel_ID_Name_Map.Name})
                warning(['On ',fn,': channel names do not match'])
            end
            FileData = varnames2struct(fileDataFields, '');
            for CHNAMESGRP = channelNamesGrouped
                chNamesGrp = CHNAMESGRP{:};
                for CHNAME = chNamesGrp
                chName = CHNAME{:};
                try
                t1 = eval([chName,'_TimeBegin']); % s
                t2 = eval([chName,'_TimeEnd']); % s
                fs = eval([chName,'_KHz'])*1000; % Hz
                t = t1:(1/fs):t2; 
                Data = eval(chName); % int
                Data = double(Data);
                if length(t) ~= length(Data)
                    err = abs(length(Data) - length(t))/length(Data);
                    warning([chName,' mismatch in time and data length by ', num2str(100*err),'%'])
                    % treat times as more accurate than sample rate 
                    t = linspace(t1, t2, length(Data));
                end

                t = seconds(t);
                GroupName = repmat(fn1, size(Data)); FileName = repmat(fn2, size(Data));
                ChannelData = varnames2struct(channelDataFields, [chName,'_']);
                Data = Data'; GroupName = GroupName'; FileName = FileName';
                T = timetable(Data,GroupName,FileName, 'RowTimes', t');
                T.Properties.Description = chName;
                T.Properties.UserData.FileData = FileData; 
                T.Properties.UserData.ChannelData = ChannelData;
                Ti = find(strcmp(channelNames, chName));
                Tbls{Ti} = [Tbls{Ti}; T];

                % if the memory is getting full, save and clear 
                sz = whos('Tbls'); sz = sz.bytes;
                if sz > sizethresh
                    save([svloc,filesep,'SaveFileH',num2str(svN),'.mat'], "Tbls","channelNames");
                    svN = svN+1;
                    clear Tbls
                    Tbls = Tbls0;
                end

                catch ME
                    warning(['On ',fn,' - ',chName,': ',ME.message])
                end
                clear chName t1 t2 t Data err GroupName FileName ChannelData T Ti
                end
            end
            clear FileData
        end
    end
end

%% run2 - consolidate saved files "vertically" 
clearvars -except channelNames svloc svN sizethresh
channelNames0 = channelNames; clear channelNames
filelist = dir(svloc);
filelist(~[filelist.isdir]);
filelist = filelist(arrayfun(@(f) f.bytes > 0, filelist));
filelist = filelist(contains({filelist.name}, 'SaveFileH'));

for chInd0 = 1:length(channelNames0)
    chName = channelNames0{chInd0};
    Tbl = [];
    svN = 1;
    for f = filelist'
        clearvars -except f filelist Tbl chName chInd0 channelNames0 ...
            svN svloc sizethresh
        fnfull = fullfile(f.folder, f.name); 
        [fp,fn,fe] = fileparts(fnfull);
        if strcmpi(fe, '.mat')
            load(fnfull);
            chInd_file = find(strcmp(chName, channelNames));
            if ~isempty(chInd_file)
                T = Tbls{chInd_file}; clear Tbls
                Tbl = [Tbl; T];

                sz = whos('Tbl'); sz = sz.bytes; 
                if sz > sizethresh
                    save([svloc,filesep,chName,'__SaveFileV',num2str(svN),'.mat'], "Tbl");
                    svN = svN+1;
                    clear Tbl
                    Tbl = [];
                end
            end
        end
    end

    save([svloc,filesep,chName,'__SaveFileV',num2str(svN),'.mat'], "Tbl");
    svN = svN+1;
    clear Tbl
end

%% helper 

function fs = getsamplerate(chName)
fs = 0;
    try 
        fs = evalinwrapper([chName,'_KHz'])*1000; % Hz
    catch ME
        warning(['On ',chName,': ',ME.message])
    end
end

function S = varnames2struct(varnames, header)
if nargin < 2
    header = '';
end
vars = cellfun(@(vn) evalinwrapper([header,vn]), varnames, 'UniformOutput', false);
S = cell2struct(vars, varnames, 2);
end

function yn = strcmpwrapper(strs1, strs2)
yn = length(strs1) == length(strs2); 
if yn 
    yn = prod(strcmp(strs1, strs2));
end
end

function var = evalinwrapper(varname)
try 
    var = evalin('base',varname);
catch ME 
    if contains(ME.message, 'Unrecognized function or variable')
        var = nan;
    else
        rethrow(ME)
    end
end
end