%% definitions 
t2rng = @(tvector) timerange(tvector(1), tvector(2));

%% user selects folder; necessary files are pulled 
folder = uigetdir; 
[~,pName] = fileparts(folder)
NSfiles = dir([folder,filesep,'*.ns*']);
NEVfiles = dir([folder,filesep,'*.nev']);
ElecXL = [dir([folder,filesep,'electrode_data.xls*']); ...
          dir([folder,filesep,pName,filesep,'electrode_data.xls*'])];

%% get electrode data 

electbl = [];
for f = ElecXL
    fTbl = readtable([f.folder,filesep,f.name]);

    % interpret coordinates
    XYZ = zeros(height(fTbl),3);
    for r = 1:height(fTbl)
        xyz = fTbl.Coordinates{r};
        xyz = sscanf(xyz,'%f');
        XYZ(r,:) = xyz';
    end
    x = XYZ(:,1); y = XYZ(:,2); z = XYZ(:,3);
    XYZ = table(x, y, z, ...
        'RowNames',fTbl.Properties.RowNames);

    fTbl = removevars(fTbl, ...
        ["Coordinates","Discard","Epileptic","OutOfBrain","Notes","LocMeeting"]);
    fTbl = [fTbl, XYZ];

    electbl = [electbl; fTbl];
end

%% get blackrock data into a time/event table 

timeBegin = datetime([inf inf inf], 'TimeZone','UTC'); 
timeEnd = datetime([-inf -inf -inf], 'TimeZone','UTC');

% NS (continuous data) files 
NStbls = {}; SampleRates = [];
for f = NSfiles
    try
    fNS = openNSx([f.folder,filesep,f.name], 'uV');
    fTbl = ns2timetable(fNS);
    timeBegin = min(timeBegin, min(fTbl.Time));
    timeEnd = max(timeEnd, max(fTbl.Time));

    % pair timetable variables with electrode info
    for c = 1:width(fTbl)
        cname = fTbl.Properties.VariableNames{c};
        if length(cname) >= 3
            cname = [cname(1),cname(3:end)]; % match naming conventions 
        end
        r = find(strcmpi(electbl.Electrode, cname));
        if ~isempty(r)
            if length(r) > 1
                warning(['Electrode name ',cname,' is not unique.']);
                r = r(1);
            end
            fTbl.Properties.VariableDescriptions{c} = electbl.Brainnetome{r};
        end
    end

    % label start/end of file 
    fTbl.Properties.Events = eventtable(...
        [min(fTbl.Time); max(fTbl.Time)], ...
        "EventLabels",...
        {['Start of File ',f.name]; ...
         ['End of File ',f.name]}, ...
        "EventEnds",[min(fTbl.Time); max(fTbl.Time)]);

    SamplingFreq = fNS.MetaTags.SamplingFreq;
    SFi = find(SampleRates == SamplingFreq);
    if isempty(SFi)
        SampleRates = [SampleRates, SamplingFreq];
        NStbls = [NStbls, {fTbl}];
    else
        try
            NStbls{SFi} = [NStbls{SFi}; fTbl];
        catch ME1
            if isequal(ME1.identifier, 'MATLAB:table:vertcat:SizeMismatch') || ...
               isequal(ME1.identifier, 'MATLAB:table:vertcat:UnequalVarNames')
                for newvarname = string(fTbl.Properties.VariableNames)
                    if ~sum(strcmpi(NStbls{SFi}.Properties.VariableNames, newvarname))
                        newvarval = nan(height(NStbls{SFi}),1);
                        eval(newvarname+" = newvarval;");
                        eval("NStbls{SFi} = addvars(NStbls{SFi},"+newvarname+");");
                    end
                end
                NStbls{SFi} = [NStbls{SFi}; fTbl];
            else
                rethrow(ME1);
            end
        end
    end
    catch ME
        warning(ME.message);
    end
end

% NEV (event data) files
NEVtbl = [];
for f = NEVfiles
    try
    fEV = openNEV([f.folder,filesep,f.name]);
    EVtbl = nev2table(fEV);
    timeBegin = min(timeBegin, min(EVtbl.Time));
    timeEnd = max(timeEnd, max(EVtbl.Time));
    NEVtbl = [NEVtbl; EVtbl];
    catch ME
        warning(ME.message);
    end
end

chnames = {};
for SFi = 1:width(NStbls)
    % pair NS with NEV
    NStbls{SFi}.Properties.Events = [NStbls{SFi}.Properties.Events; NEVtbl];
    % get all channel names 
    chnames = [chnames, NStbls{SFi}.Properties.VariableNames];
end
chnames = unique(chnames); chnames = string(chnames);

%% user selects channels 
channelIndexRec = listdlg("PromptString","Recording Channel(s)", ...
    "ListString",chnames, "SelectionMode","multiple"); 
channelIndexStim = listdlg("PromptString","Stimulus Channel(s)", ...
    "ListString",chnames, "SelectionMode","multiple"); 
channelIndexTrig = listdlg("PromptString","Stimulus TRIGGER Channel(s)", ...
    "ListString",chnames, "SelectionMode","multiple"); 
channelIndexInspect = listdlg("PromptString","Inspect Channel(s)", ...
    "ListString",chnames, "SelectionMode","multiple"); 
for chindlist = ["Rec", "Stim", "Trig", "Inspect"]
    eval("channelName"+chindlist+" = chnames(channelIndex"+chindlist+");");
end

%% detect triggers 

% remove trigger pulse trains and convert to event table 
trigtbl = [];
for chname = channelNameTrig
    lbl = "Trigger "+chname;
    for SFi = 1:width(NStbls)
        NStbl = NStbls{SFi};
        if sum(strcmp(NStbl.Properties.VariableNames, chname))
            t = NStbl.Time; trig = NStbl.(chname);
            trig = [false; diff(trig) > 1e3]; % rising edge
            t = t(trig);
            trigtbl = [trigtbl; eventtable(t, "EventLabels",lbl, "EventEnds",t)];
            NStbl = removevars(NStbl, chname);
            NStbls{SFi} = NStbl;
        end
    end
end

% group triggers by time 
trigboundtime = 1; % sec
trigtime = trigtbl.Time; % all trigger types!
trigbound = seconds(diff(trigtime)) > trigboundtime;
trigstart = trigtime([trigbound; false]); 
trigend = trigtime([false; trigbound]);
trigstart = [trigtime(1); trigstart];
trigend = [trigend; trigtime(end)];
trigwinds = [trigstart, trigend];

for SFi = 1:width(NStbls)
    NStbls{SFi}.Properties.Events = [NStbls{SFi}.Properties.Events; trigtbl];
end

%% construct main table 
% everything resampled to the same time rows
%if length(SampleRates) > 1
[~,SFi] = min(abs(SampleRates-1000));
SamplingFreq = SampleRates(SFi);
SamplingFreq = inputdlg('Select main frequency to resample:', ...
    'main table sampling selection', 1, {num2str(SamplingFreq)});
SamplingFreq = eval(SamplingFreq{1});
[~,SFi] = min(abs(SampleRates-SamplingFreq));
MainTable = NStbls{SFi};
MainTable = myRetime(MainTable, SampleRates(SFi), nan);
MainTable = retime(MainTable, 'regular', 'nearest', 'SampleRate',SamplingFreq);
SFj = true(size(SampleRates)); SFj(SFi) = false; SFj = find(SFj);
for SFi = SFj
    tbl = myRetime(NStbls{SFi}, SampleRates(SFi), nan);
    SampleRatio = ceil(SampleRates(SFi)/SamplingFreq);
    if SampleRatio > 1
        tbl.Variables = movmean(tbl.Variables, SampleRatio, 'omitnan');
    end
    tbl = retime(tbl, MainTable.Time, 'nearest');
    MainTable = [MainTable, tbl];
end
%end

%% determine the baseline 
% can any of this be changed to use MainTable? 

% assess file start/end, packet loss, etc
discont = [];
for SFi = 1:width(NStbls)
    evtbl = NStbls{SFi}.Properties.Events;
    for disctype = ["Start of file", "End of file", "SerialDigitalIO"]
        discevt = contains(evtbl.EventLabels, disctype);
        discont = [discont; evtbl.Time(discevt)];
    end
    disctype = "Data Loss";
        discevt = contains(evtbl.EventLabels, disctype);
        discont = [discont; evtbl.EventEnds(discevt)];
    inan = isnan(NStbls{SFi}.Variables); inan = sum(inan,2);
    discont = [discont; NStbls{SFi}.Time(find(inan))];
end

% canditate periods in between triggers, etc
candwinds = [];
boundtime = seconds(.4); % sec 
trigend_ = [timeBegin; trigend];
trigstart_ = [trigstart; timeEnd];
for itrig = 1:length(trigstart_)
    t1 = trigend_(itrig)+boundtime; t2 = trigstart_(itrig)-boundtime;
    disci = (discont > t1) & (discont < t2);
    disci = discont(disci);
    if ~isempty(disci)        
        for d = disci'
            candwinds = [candwinds; t1, d-boundtime];
            t1 = d+boundtime;
        end
    end
    candwinds = [candwinds; t1, t2];
end
candwindslong = (candwinds(:,2) - candwinds(:,1)) >= minutes(3); % MIN DURATION
candwinds = candwinds(candwindslong,:);

% outlier detection 
OLtbls = cell(size(NStbls)); 
for SFi = 1:width(NStbls)
    io = isoutlier(NStbls{SFi}, 'mean');
    numOL = sum(io,2);
    OLtbls{SFi} = timetable(NStbls{SFi}.Time, numOL);
end

% outliers during candidate windows 
OLwinds = zeros(height(candwinds),1);
for iwind = 1:height(OLwinds)
    trng = timerange(candwinds(iwind,1), candwinds(iwind,2));
    t1 = timeEnd; t2 = timeBegin;
    for SFi = 1:width(NStbls)
        io = OLtbls{SFi}(trng,:);
        t1 = min(t1, min(io.Time)); t2 = max(t2, max(io.Time));
        OLwinds(iwind) = OLwinds(iwind) + sum(io.Variables.^2);
    end
    OLwinds(iwind) = OLwinds(iwind)/seconds(t2-t1);
end

%% user confirms baseline 
[~,iwind] = min(OLwinds);
fig1 = figure; [~,hAXs] = myStackedPlot(MainTable, ...
    [channelNameRec, channelNameStim, channelNameInspect]);
trngBaseline = inputdlg({'Start', 'End'}, 'BASELINE Segment', ...
    1, string(candwinds(iwind,:)));
trngBaseline = datetime(trngBaseline, 'TimeZone',candwinds.TimeZone);
tblBaselineMain = mySelect(MainTable, t2rng(trngBaseline));
tblsBaseline = cellfun(@(T) mySelect(T,t2rng(trngBaseline)), ...
    NStbls, 'UniformOutput',false);
tblsBaseline = [{tblBaselineMain}, tblsBaseline];

% nan check
[trngBaseline,tblsBaseline] = nancheck(trngBaseline,tblsBaseline,hAXs);
tblBaselineMain = tblsBaseline{1};

% indicate on plot
xrng = trngBaseline;
xrng = [xrng, xrng]'; xrng = xrng(:)';
xrng = [xrng, xrng(1)]; 
for hAX = hAXs'
    yrng = ylim(hAX); 
    yrng = [yrng, fliplr(yrng)];
    yrng = [yrng, yrng(1)];
    patch(hAX, xrng, yrng, [0,0,1], 'FaceAlpha',.5);
end

%% user confirms marked stimulation 

%% show channel stats to reject noisy channels 

%% helpers 

function Tbl = mySelect(Tbl, times)
Tbl = Tbl(times,:);
Tbl.Properties.Events = Tbl.Properties.Events(times,:);
end

function [trng, tbls] = nancheck(trng, tbls, hAXs)
tnan = [];
for SFi = 1:width(tbls)
    TBLi = tbls{SFi};
    inan = sum(isnan(TBLi.Variables),2);
    tnan = [tnan; TBLi.Time(find(inan))];
end
if ~isempty(tnan)
    for hAX = hAXs'
        plot(hAX, tnan, zeros(size(tnan)), '*r');
    end
    xlim(trng)
    trimsel = questdlg('Selection contains missing/nan values.', ...
        'Reselect due to nan values', ...
        'Start after nan', 'End before nan', 'Proceed anyway', ...
        'Proceed anyway');
    if isempty(trimsel)
        error('Selection must be made.')
    elseif strcmp(trimsel, 'Start after nan')
        trngBasleline(1) = max(tnan) + seconds(.001);
    elseif strcmp(trimsel, 'End before nan')
        trngBasleline(2) = min(tnan) - seconds(.001);
    end
    if ~strcmp(trimsel, 'Proceed anyway')
        tbls = cellfun(@(T) mySelect(T,t2rng(trngBasleline)), ...
            tbls, 'UniformOutput',false);
    end
end
end