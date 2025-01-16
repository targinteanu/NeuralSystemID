function TT = ns2timetable(nsOrFilename)

if nargin < 1
    nsOrFilename = [];
end
cl = class(nsOrFilename);

%% load file

if isstruct(nsOrFilename)
    NS = nsOrFilename;
else
    if isempty(nsOrFilename)
        [fn,fp] = uigetfile({'*.ns*'; '*.mat'});
        [~,fn,fe] = fileparts(fn);
    elseif ischar(nsOrFilename) || isstring(nsOrFilename)
        [fp,fn,fe] = fileparts(nsOrFilename);
    end
    if strcmpi(fe,'.mat')
        load(fullfile(fp,[fn,fe]), 'NS2', 'ns2', 'NS5', 'ns5', 'NS', 'NS');
        for vtry = {'NS2', 'ns2', 'NS5', 'ns5', 'NS'}
            if exist(vtry{:}, 'var')
                NS = eval(vtry{:});
            end
        end
    else
        openNSx(fullfile(fp,[fn,fe]), 'uV');
        NS = eval(['NS',fe(end)]);
    end
    clear vtry
end

%% interpret data from loaded file 
% obtain usable data variables and other information from the file. 

% Get timing data:
SamplingFreq = NS.MetaTags.SamplingFreq;
datalen = NS.MetaTags.DataPoints;
try
    t1Rel = NS.MetaTags.Timestamp / NS.MetaTags.TimeRes;
catch
    t1Rel = 0;
end
if isfield(NS.MetaTags, 'DataDurationSec')
    datadur = NS.MetaTags.DataDurationSec; 
elseif isfield(NS.MetaTags, 'DataPointsSec')
    datadur = NS.MetaTags.DataPointsSec;
else
    datadur = datalen/SamplingFreq;
end

% handle packet loss 
if (length(t1Rel) > 1) || (length(datadur) > 1)
    warning('Packet loss. Largest packet will be selected.') % replace with something better?
    tStartEnd = [t1Rel; t1Rel + datadur];
    packetdur = diff(tStartEnd);
    [~,p] = max(packetdur); 
    dta = NS.Data{p}; datalen = datalen(p); datadur = datadur(p);
    if length(t1Rel) > 1
        t1Rel = t1Rel(p);
    end
else
    dta = NS.Data;
end

% tRel = time relative to start of recording, in seconds 
% t = absolute date/time 
% t0 = date/time at start of recording 
tRel = linspace(0, datadur, datalen) + t1Rel;
t = seconds(tRel);
t0 = datetime(NS.MetaTags.DateTime); 
t = t+t0; 

lbl = {NS.ElectrodesInfo.Label};
lbl = upper(lbl);
unitname = {NS.ElectrodesInfo.AnalogUnits};

%% construct table

X = dta'; X = double(X);

TT = array2timetable(X,"RowTimes",t,"VariableNames",lbl); 
TT.Properties.VariableUnits = unitname;
%TT.Properties.SampleRate = NS.MetaTags.SamplingFreq;

TT.Properties.Description = [...
    'file ',NS.MetaTags.FilePath,filesep,NS.MetaTags.Filename,NS.MetaTags.FileExt,...
    ', comment ',NS.MetaTags.Comment];

% sample rate check 
fs_Signal = NS.MetaTags.SamplingFreq;
dT = seconds(diff(t));
dTmax = max(dT); dTmin = min(dT);
fsmax = 1/dTmin; fsmin = 1/dTmax;
errthresh = .01;
if abs((fsmax-fs_Signal)/fs_Signal) > errthresh
    warning('Sample Rate is incorrectly reported or inconsistent.')
end
if abs((fsmin-fs_Signal)/fs_Signal) > errthresh
    warning('Sample Rate is incorrectly reported or inconsistent.')
end

end