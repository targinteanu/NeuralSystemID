%% setup
filelist = dir(uigetdir);
file1 = filelist(~[filelist.isdir]);
file1 = file1(arrayfun(@(f) f.bytes > 0, file1));
file1 = file1(1); load(fullfile(file1.folder, file1.name));
%{
channelNames = {...
    'ANALOG_IN_1', ...
    'ANALOG_IN_2', ...
    'LFP_01', ...
    'LFP_02', ...
    'Macro_LFP_01', ...
    'Macro_LFP_02', ...
    'Macro_RAW_01', ...
    'Macro_RAW_02', ...
    'RAW_01', ...
    'RAW_02', ...
    'SEG_01', ...
    'SEG_02', ...
    'SPK_01', ...
    'SPK_02' ...
    };
%}
channelNames = {Channel_ID_Name_Map.Name};
Tbls = cell(size(channelNames));
fileDataFields = {...
    'SF_DRIVE_CONF', ...
    ... 'SF_DRIVE_DOWN', ...
    'SF_DRIVE_SET_POS', ...
    ... 'SF_DRIVE_STOP', ...
    ... 'SF_DRIVE_UP', ...
    'SF_LEVEL' ...
    }; 
channelDataFields = {'BitResolution', 'Gain'};

%% run
for f = filelist'
    clearvars -except Tbls channelNames channelDataFields fileDataFields f filelist
    if (~f.isdir) && (f.bytes > 0)
        fnfull = fullfile(f.folder, f.name); 
        [fp,fn,fe] = fileparts(fnfull);
        fni = find(fn == '.');
        fn1 = fn(1:(fni-1)); fn2 = fn((fni+1):end);
        fn1 = string(fn1); fn2 = string(fn2);
        if strcmpi(fe, '.mat')
            load(fnfull)
            FileData = varnames2struct(fileDataFields, '');
            for ID_name = Channel_ID_Name_Map'
                chName = ID_name.Name;
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
                T = timetable(Data',GroupName',FileName', 'RowTimes', t');
                T.Properties.Description = chName;
                T.Properties.UserData.FileData = FileData; 
                T.Properties.UserData.ChannelData = ChannelData;
                Ti = find(strcmp(channelNames, chName));
                Tbls{Ti} = [Tbls{Ti}; T];
                catch ME
                    warning(['On ',fn,' - ',chName,': ',ME.message])
                end
                clear chName t1 t2 t Data err GroupName FileName ChannelData T Ti
            end
            clear FileData
        end
    end
end

%% helper 
function S = varnames2struct(varnames, header)
if nargin < 2
    header = '';
end
vars = cellfun(@(vn) evalin('base',[header,vn]), varnames, 'UniformOutput', false);
S = cell2struct(vars, varnames, 2);
end