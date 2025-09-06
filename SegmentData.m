%% user selects folder; necessary files are pulled 
disp('Please select subject folder...')
folder = uigetdir; 
[~,pName] = fileparts(folder)
NSfiles = dir([folder,filesep,'*.ns*']);
NEVfiles = dir([folder,filesep,'*.nev']);
ElecXL = [dir([folder,filesep,'electrode_data.xls*']); ...
          dir([folder,filesep,pName,filesep,'electrode_data.xls*'])];
thisfilename = mfilename("fullpath");

%% get electrode data 
disp('Reading electrode localization data...')
if isempty(ElecXL)
    warning('No electrode localization data found!')
end

electbl = [];
for f = ElecXL'
    fTbl = readtable([f.folder,filesep,f.name]);

    % interpret coordinates
    if sum(strcmpi(fTbl.Properties.VariableNames, 'coordinates'))
    XYZ = zeros(height(fTbl),3);
    for r = 1:height(fTbl)
        xyz = fTbl.Coordinates{r};
        xyz = sscanf(xyz,'%f');
        XYZ(r,:) = xyz';
    end
    x = XYZ(:,1); y = XYZ(:,2); z = XYZ(:,3);
    XYZ = table(x, y, z, ...
        'RowNames',fTbl.Properties.RowNames);
    else
        XYZ = [];
    end

    fTbl = removevarswrapper(fTbl, ...
        ["Coordinates","Discard","Epileptic","OutOfBrain","Notes","LocMeeting"]);
    fTbl = [fTbl, XYZ];

    electbl = [electbl; fTbl];
end

clear x y z XYZ xyz f fTbl

%% get blackrock data into a time/event table 
disp('Accessing blackrock data...')

timeBegin = datetime([inf inf inf], 'TimeZone','UTC'); 
timeEnd = datetime([-inf -inf -inf], 'TimeZone','UTC');

% NS (continuous data) files 
disp('  - continuous data...')
tbls = {}; SampleRates = [];
for f = NSfiles'
    try
    fNS = openNSx([f.folder,filesep,f.name], 'uV');
    fTbl = ns2timetable(fNS);

    disp(['File: ',f.name,...
          ' starts at ',char(min(fTbl.Time)),' and',...
          ' ends at ',char(max(fTbl.Time))]);
    varnames = fTbl.Properties.VariableNames;
    varnames = cellfun(@(v) [v,'|'], varnames, 'UniformOutput', false);
    disp([varnames{:}]);
    proc = input('Include this file? (Y/N) ', "s");
    if strcmpi(proc, 'Y')

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
        tbls = [tbls, {fTbl}];
    else
        try
            tbls{SFi} = [tbls{SFi}; fTbl];
        catch ME1
            if isequal(ME1.identifier, 'MATLAB:table:vertcat:SizeMismatch') || ...
               isequal(ME1.identifier, 'MATLAB:table:vertcat:UnequalVarNames')
                for newvarname = string(fTbl.Properties.VariableNames)
                    if ~sum(strcmpi(tbls{SFi}.Properties.VariableNames, newvarname))
                        newvarval = nan(height(tbls{SFi}),1);
                        eval(newvarname+" = newvarval;");
                        eval("tbls{SFi} = addvars(tbls{SFi},"+newvarname+");");
                    end
                end
                for newvarname = string(tbls{SFi}.Properties.VariableNames)
                    if ~sum(strcmpi(fTbl.Properties.VariableNames, newvarname))
                        newvarval = nan(height(fTbl),1);
                        eval(newvarname+" = newvarval;");
                        eval("fTbl = addvars(fTbl,"+newvarname+");");
                    end
                end
                tbls{SFi} = [tbls{SFi}; fTbl];
            else
                rethrow(ME1);
            end
        end
    end
    end
    catch ME
        warning(ME.message);
    end
end

% NEV (event data) files
disp('  - event data...')
NEVtbl = [];
for f = NEVfiles'
    try
    fEV = openNEV([f.folder,filesep,f.name], 'nosave');
    EVtbl = nev2table(fEV);

    disp(['File: ',f.name,...
          ' first event at ',char(min(EVtbl.Time)),' and',...
          ' last event at ',char(max(EVtbl.Time))]);
    disp(unique(EVtbl.EventLabels));
    proc = input('Include this file? (Y/N) ', "s");
    if strcmpi(proc, 'Y')

    timeBegin = min(timeBegin, min(EVtbl.Time));
    timeEnd = max(timeEnd, max(EVtbl.Time));
    NEVtbl = [NEVtbl; EVtbl];

    end
    catch ME
        warning(ME.message);
    end
end

% ensure time ascending 
if numel(NEVtbl)
    [~,sortind] = sort(NEVtbl.Time);
    NEVtbl = NEVtbl(sortind,:);
end

chnames = {};
for SFi = 1:width(tbls)
    if numel(tbls{SFi})
    % ensure time ascending 
    [~,sortind] = sort(tbls{SFi}.Time);
    tbls{SFi} = tbls{SFi}(sortind,:);
    % pair NS with NEV
    tbls{SFi}.Properties.Events = [tbls{SFi}.Properties.Events; NEVtbl];
    % get all channel names 
    chnames = [chnames, tbls{SFi}.Properties.VariableNames];
    end
end
chnames = unique(chnames); chnames = string(chnames);

clear f fNS fEV fTbl EVtbl c r cname proc varnames sortind

%% user selects channels 
disp('Please specify channels...')
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
clear chindlist

%% detect triggers 
disp('Detecting triggers...')

% ensure trigger signal is not duplicated 
for chname = channelNameTrig
    namefound = false;
    for SFi = 1:width(tbls)
        namefound_ = sum(strcmp(chname, tbls{SFi}.Properties.VariableNames));
        if namefound && namefound_
            error(['Trigger signal ',char(chname),...
                ' is duplicated in multiple recording files!']);
        else
            namefound = namefound || namefound_;
        end
    end
end

% remove trigger pulse trains and convert to event table 
trigtbl = [];
for chname = channelNameTrig
    lbl = "Trigger "+chname;
    for SFi = 1:width(tbls)
        NStbl = tbls{SFi};
        if sum(strcmp(NStbl.Properties.VariableNames, chname))
            t = NStbl.Time; trig = NStbl.(chname);
            trig = [false; diff(trig) > 1e3]; % rising edge
            t = t(trig);
            trigtbl = [trigtbl; eventtable(t, "EventLabels",lbl, "EventEnds",t)];
            NStbl = removevars(NStbl, chname);
            if isempty(NStbl)
                SampleRates(SFi) = nan;
            else
                tbls{SFi} = NStbl;
            end
        end
    end
end
tbls = tbls(~isnan(SampleRates)); SampleRates = SampleRates(~isnan(SampleRates));

% group triggers by time 
trigboundtime = 10; % sec (defines a burst of triggers)
trigtime = trigtbl.Time; % all trigger types!
trigwinds = groupbytime(trigtime, trigboundtime, timeBegin);
trigstart = trigwinds(:,1); trigend = trigwinds(:,2);

% describe triggers 
trignames = repmat("", height(trigwinds),1);
for itrig = 1:height(trigwinds)
    trigname = trignames(itrig);
    trigtbl_ = trigtbl(timerange(trigwinds(itrig,1),trigwinds(itrig,2)),:);
    trignames_ = unique(trigtbl_.EventLabels);
    for trigname_ = trignames_'
        trigname = trigname+"; "+trigname_;
    end
    trigname = char(trigname);
    trignames(itrig) = string(trigname(3:end));
end

for SFi = 1:width(tbls)
    tbls{SFi}.Properties.Events = [tbls{SFi}.Properties.Events; trigtbl];
end

clear NStbl t trig lbl trigname trigtbl_ trignames_ trigname_ itrig
clear namefound namefound_

%% construct main table 
% everything resampled to the same time rows
disp('Resampling...')
%if length(SampleRates) > 1
[~,SFi] = min(abs(SampleRates-1000));
SamplingFreq = SampleRates(SFi);
SamplingFreq = inputdlg('Select main frequency to resample:', ...
    'main table sampling selection', 1, {num2str(SamplingFreq)});
SamplingFreq = eval(SamplingFreq{1});
[~,SFi] = min(abs(SampleRates-SamplingFreq));
MainTable = tbls{SFi};
MainTable = myRetime(MainTable, SampleRates(SFi), nan);
MainTable = retime(MainTable, 'regular', 'nearest', 'SampleRate',SamplingFreq);
SFj = true(size(SampleRates)); SFj(SFi) = false; SFj = find(SFj);
for SFi = SFj
    tbl = myRetime(tbls{SFi}, SampleRates(SFi), nan);
    SampleRatio = ceil(SampleRates(SFi)/SamplingFreq);
    if SampleRatio > 1
        tbl.Variables = movmean(tbl.Variables, SampleRatio, 'omitnan');
    end
    tbl = retime(tbl, MainTable.Time, 'nearest');
    MainTable = [MainTable, tbl];
end
%end

chnames = unique(MainTable.Properties.VariableNames);

clear chname

%% determine the baseline 
% can any of this be changed to use MainTable? 
disp('Detecting baseline...')

% assess file start/end, packet loss, etc
discont = [];
for SFi = 1:width(tbls)
    evtbl = tbls{SFi}.Properties.Events;
    for disctype = ["Start of File", "End of File", "SerialDigitalIO"]
        discevt = contains(evtbl.EventLabels, disctype);
        discont = [discont; evtbl.Time(discevt)];
    end
    disctype = "Data Loss";
        discevt = contains(evtbl.EventLabels, disctype);
        discont = [discont; evtbl.EventEnds(discevt)];
    inan = isnan(tbls{SFi}.Variables); inan = sum(inan,2);
    inantime = tbls{SFi}.Time(find(inan));
    timetol = 3/SampleRates(SFi); % dur to be considered same event (s)
    inantime = inantime([true; seconds(diff(inantime)) > timetol]);
    discont = [discont; inantime];
end

clear disctype discevt inan evtbl timetol inantime

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

clear trigstart_ trigend_ t1 t2 disci d discont

% outlier detection 
OLtbls = cell(size(tbls)); 
for SFi = 1:width(tbls)
    io = isoutlier(tbls{SFi}, 'mean');
    numOL = sum(io,2);
    OLtbls{SFi} = timetable(tbls{SFi}.Time, numOL);
end

% outliers during candidate windows 
OLwinds = zeros(height(candwinds),1);
for iwind = 1:height(OLwinds)
    trng = timerange(candwinds(iwind,1), candwinds(iwind,2));
    t1 = timeEnd; t2 = timeBegin;
    for SFi = 1:width(tbls)
        io = OLtbls{SFi}(trng,:);
        t1 = min(t1, min(io.Time)); t2 = max(t2, max(io.Time));
        OLwinds(iwind) = OLwinds(iwind) + sum(io.Variables.^2);
    end
    OLwinds(iwind) = OLwinds(iwind)/seconds(t2-t1);
end

clear t1 t2 numOL io trng iwind

%% user confirms baseline 
disp('Please confirm baseline condition...')
[~,iwind] = min(OLwinds);
fig1 = figure('Units','normalized', 'Position',[.05,.05,.9,.9]); 
[~,hAXs] = myStackedPlot(MainTable, ...
    [channelNameRec, channelNameStim, channelNameInspect]);
trngBaseline = candwinds(iwind,:);

doAgain = true;
while doAgain
trngBaseline = inputdlg({'Start', 'End'}, 'BASELINE Segment', ...
    1, string(trngBaseline));
trngBaseline = datetime(trngBaseline, 'TimeZone',candwinds.TimeZone);
tblBaselineMain = mySelect(MainTable, trngBaseline, true);
tblsBaseline = cellfun(@(T) mySelect(T, trngBaseline, true), ...
    tbls, 'UniformOutput',false);
tblsBaseline = [{tblBaselineMain}, tblsBaseline];

% nan check
[trngBaseline,tblsBaseline,channelNameReject] = ...
    nancheck(trngBaseline,tblsBaseline,hAXs);
if numel(channelNameReject)
    for chname = channelNameReject
        chnames = chnames(~strcmp(chnames, chname));
    end
end
clear chname
tblBaselineMain = tblsBaseline{1};

% indicate on plot
hPat = plottrng(trngBaseline, [0,0,1], hAXs);
pause(.01); drawnow; pause(.01);

doAgainQst = questdlg('Accept?', 'Accept Time Range', ...
    'Yes', 'Try Again', 'Stop Program', 'Yes');
if strcmp(doAgainQst, 'Yes')
    doAgain = false;
elseif strcmp(doAgainQst, 'Try Again')
    delete(hPat);
else
    error('Terminated by user.')
end
end

clear doAgain doAgainQst hPat

%% show channel stats
% altered version of tblChannelSummary.m
disp('Calculating baseline characteristics...')
BaselineData = tblBaselineMain.Variables;
varnames = tblBaselineMain.Properties.VariableNames; 
vardescs = tblBaselineMain.Properties.VariableDescriptions;
for d = 1:length(vardescs)
    if length(vardescs{d}) > 20
        vardescs{d} = [vardescs{d}(1:3),'...',vardescs{d}((end-20):end)];
    end
end

% ignore nan
BaselineData_inan = sum(isnan(BaselineData),2);
BaselineData = BaselineData(~BaselineData_inan, :);

% band(s) power and range 
disp('  - band power:')
BaselineSD = std(BaselineData);
BaselineBetaPower = bandpower(BaselineData, SamplingFreq, [13,30]);
BaselineThetaPower = bandpower(BaselineData, SamplingFreq, [4,9]);
BaselineGammaPower = bandpower(BaselineData, SamplingFreq, [30,80]);

% # of outliers 
disp('  - outlier detection:')
io = isoutlier(BaselineData, 'mean');
BaselineNumOL = sum(io,1);

% display format selection 
toshowas = [false, false]; % stem, grid
showas = questdlg('Show channel data as:', 'Select data display format', ...
    'stem', 'grid', 'both', 'both');
if strcmp(showas, 'stem')
    toshowas(1) = true;
end
if strcmp(showas, 'grid')
    toshowas(2) = true;
end
if strcmp(showas, 'both')
    toshowas(:) = true;
end

% stem display
if toshowas(1)
    fig2 = figure('Units','normalized', 'Position',[.05,.05,.9,.9]);
    subplot(6,1,1);
    plotstem(BaselineThetaPower, varnames); 
    ylabel('\theta Power'); grid on; axis tight;
    xl = xlim();
    subplot(6,1,2);
    plotstem(BaselineBetaPower,  varnames); 
    ylabel('\beta Power'); grid on; axis tight;
    subplot(6,1,3);
    plotstem(BaselineGammaPower, varnames); 
    ylabel('\gamma Power'); grid on; axis tight;
    subplot(6,1,4);
    plotstem(BaselineSD, varnames); 
    ylabel('S.D.'); grid on; axis tight;
    subplot(6,1,5);
    plotstem(BaselineNumOL, varnames); 
    ylabel('# OutL.'); grid on; axis tight;
    subplot(6,1,6);
    text(1:length(varnames), zeros(size(varnames)), vardescs, ...
        'Rotation',90, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
    xlim(xl); xticks([]); yticks([]);
    sgtitle('Baseline Properties')
end

% grid display 
if toshowas(2)
    fig3 = figure('Units','normalized', 'Position',[.05,.05,.9,.9]);
    subplot(1,5,1); 
    plotgrid(BaselineThetaPower, varnames, 21, 3); 
    title('\theta Power'); axis tight;
    subplot(1,5,2); 
    plotgrid(BaselineBetaPower,  varnames, 21, 3); 
    title('\beta Power'); axis tight;
    subplot(1,5,3); 
    plotgrid(BaselineGammaPower, varnames, 21, 3); 
    title('\gamma Power'); axis tight;
    subplot(1,5,4); 
    plotgrid(BaselineSD, varnames, 21, 3); 
    title('S.D.'); axis tight;
    subplot(1,5,5); 
    plotgrid(BaselineNumOL, varnames, 21, 3); 
    title('# OutL.'); axis tight;
    sgtitle('Baseline Properties')
end

pause(.01); drawnow; pause(.01);

clear showas xl

%% inspect and reject noisy channels 
disp('Please inspect and specify noisy channels to reject...')

% final inspection 
channelIndexInspect2 = listdlg("PromptString","Inspect Channel(s)", ...
    "ListString",chnames, "SelectionMode","multiple"); 
channelNameInspect2 = chnames(channelIndexInspect2);
if sum(channelIndexInspect2)
    fig4 = figure('Units','normalized', 'Position',[.05,.05,.9,.9]); 
    myStackedPlot(MainTable, channelNameInspect2);
end

% rejection 
channelIndexReject = listdlg("PromptString","REJECT Channel(s)", ...
    "ListString",chnames, "SelectionMode","multiple"); 
channelNameReject = [channelNameReject, chnames(channelIndexReject)];
MainTable = removevars(MainTable, channelNameReject);
tblBaselineMain = removevars(tblBaselineMain, channelNameReject);
tbls = cellfun(@(tbl) myRemoveVars(tbl, channelNameReject), ...
    tbls, 'UniformOutput',false);
tblsBaseline = cellfun(@(tbl) myRemoveVars(tbl, channelNameReject), ...
    tblsBaseline, 'UniformOutput',false);

%% user confirms marked stimulation 
disp('Please confirm stimulation triggers...')
figure(fig1);
trngTrig = trigwinds'; comTrig = trignames;
[trngTrig, comTrig, tblTrigMain, tblsTrig] = ...
    usrconfirm(trngTrig, comTrig, MainTable, tbls, hAXs, [1,0,0]);

%% determine unmarked stimulation 
disp('Please confirm stimulation without triggers...')

% initial screen for outliers 
stimpoints = detectStimPoints(MainTable.Variables', tblBaselineMain.Variables', .1);
stimtime = MainTable.Time(stimpoints);

% exclude baseline and marked stim 
stimtime = stimtime( (stimtime<min(trngBaseline)) | (stimtime>max(trngBaseline)));
for itrig = 1:width(trngTrig)
    stimtime = stimtime( ...
        (stimtime<min(trngTrig(:,itrig))) | ...
        (stimtime>max(trngTrig(:,itrig))));
end

% append appropriate events 
stimevt = eventtable(stimtime, "EventEnds",stimtime, ...
    "EventLabels","Presumed Stim");
MainTable.Properties.Events = [MainTable.Properties.Events; stimevt];
for SFi = 1:width(tbls)
    tbls{SFi}.Properties.Events = [tbls{SFi}.Properties.Events; stimevt];
end

% group stim by time 
stimboundtime = 1; % sec (defines a burst of triggers)
stimwinds = groupbytime(stimtime, stimboundtime, timeBegin);

% user confirms 
figure(fig1);
trngStimNoTrig = stimwinds'; comStim = repmat("Presumed Stim", height(stimwinds), 1);
[trngStimNoTrig, comStim, tblStimNoTrigMain, tblsStimNoTrig] = ...
    usrconfirm(trngStimNoTrig, comStim, MainTable, tbls, hAXs, [1,0,1]);

clear itrig

%% segment serial comm sessions 
disp('Please confirm serial events...')

% search for serial comm 
evtbl = MainTable.Properties.Events;
srlevt = contains(evtbl.EventLabels, "SerialDigitalIO");
srltime = evtbl.Time(srlevt);

% separate out marked stimulation 
srltimeMarkedStim = []; 
for itrig = 1:width(trngTrig)
    stimtime = ...
        (srltime > min(trngTrig(:,itrig))) & ...
        (srltime < max(trngTrig(:,itrig)));
    srltimeMarkedStim = [srltimeMarkedStim; srltime(stimtime)];
    srltime = srltime(~stimtime);
end
% separate out unmarked stimulation 
srltimeUnmarkedStim = []; 
for istim = 1:width(trngStimNoTrig)
    stimtime = ...
        (srltime > min(trngStimNoTrig(:,istim))) & ...
        (srltime < max(trngStimNoTrig(:,istim)));
    srltimeUnmarkedStim = [srltimeUnmarkedStim; srltime(stimtime)];
    srltime = srltime(~stimtime);
end

% group by time 
srlboundtime = 15; % s (defines a session of serial comms) 
% no stim
srlwinds = groupbytime(srltime, srlboundtime, timeBegin);
% marked stim 
srlwindsMarkedStim = groupbytime(srltimeMarkedStim, srlboundtime, timeBegin);
% unmarked stim 
srlwindsUnmarkedStim = groupbytime(srltimeUnmarkedStim, srlboundtime, timeBegin);
% recombine 
% alternatively, should these all be processed and saved separately?
srlwinds = [srlwinds; srlwindsMarkedStim; srlwindsUnmarkedStim];

% describe srl 
srlnames = repmat("", height(srlwinds),1);
for isrl = 1:height(srlwinds)
    srlname = srlnames(isrl);
    srltbl_ = evtbl(timerange(srlwinds(isrl,1),srlwinds(isrl,2)),:);
    srlnames_ = unique(srltbl_.EventLabels);
    for srlname_ = srlnames_'
        srlname = srlname+"; "+srlname_;
    end
    srlname = char(srlname);
    srlnames(isrl) = string(srlname(3:end));
end
for isrl = 1:height(srlnames)
    if isrl > height(srlnames) - height(srlwindsUnmarkedStim)
        srlnames(isrl) = "Presumed Stim; "+srlnames(isrl);
    elseif isrl > height(srlnames) - height(srlwindsUnmarkedStim) - height(srlwindsMarkedStim)
        srlnames(isrl) = "Triggered Stim; "+srlnames(isrl);
    else
        srlnames(isrl) = "No Stim; "+srlnames(isrl);
    end
end

% user confirms 
figure(fig1);
trngSrl = srlwinds'; comSrl = srlnames;
[trngSrl, comSrl, tblSrlMain, tblsSrl] = ...
    usrconfirm(trngSrl, comSrl, MainTable, tbls, hAXs, [0,1,0]);

clear srlname srltbl_ srlnames_ srlname_ isrl istim

%% segment misc 
disp('Please label any other data to save...')
% user confirms 
figure(fig1);
trngMisc = repmat(datetime(NaT,'TimeZone',timeBegin.TimeZone), 1,2)';
comMisc = "Segment any miscellaneous time periods, or close this window to cancel.";
[trngMisc, comMisc, tblMiscMain, tblsMisc] = ...
    usrconfirm(trngMisc, comMisc, MainTable, tbls, hAXs, [1,1,0]);

%% save all
disp('Saving segmented data...')

thisfilever = getFileVersion(thisfilename);
[~,thisfilename] = fileparts(thisfilename);
svname = [pName,'_',thisfilename,'_',thisfilever];

tbls = [{MainTable}, tbls];

% make save folder 
folder = [folder,filesep,svname];
mkdir(folder);
svname = fullfile(folder, svname);

% save vars (table lists)
disp('  - data file:')
save(svname, ...
    ..."tbls", ...
    "SampleRates", ...
    "tblsMisc", "tblsSrl", "tblsStimNoTrig", "tblsTrig", "tblsBaseline", ...
    "channelNameRec", "channelNameStim", "channelNameTrig", "channelNameReject", ...
    "trngMisc", "trngSrl", "trngStimNoTrig", "trngTrig", "trngBaseline", ...
    "-v7.3");

% save fig1 and fig4 
disp('  - figure 1:')
saveasmultiple(fig1, [svname,'_TimePlot1']);
if sum(channelIndexInspect2)
    disp('  - figure 4:')
    saveasmultiple(fig4, [svname,'_TimePlot2']);
end
if toshowas(1)
    % save fig2
    disp('  - figure 2:')
    saveasmultiple(fig2, [svname,'_StemPlot']);
end
if toshowas(2)
    % save fig3
    disp('  - figure 3:')
    saveasmultiple(fig3, [svname,'_GridPlot']);
end

disp('...Done!')

%% helpers 


function Tbl = removevarswrapper(Tbl, varslist)
varspresent = true(size(varslist));
for v = 1:length(varslist)
    var = varslist(v);
    if ~sum(strcmp(Tbl.Properties.VariableNames, var))
        warning("Unrecognized table variable name: "+var);
        varspresent(v) = false;
    end
end
varslist = varslist(varspresent);
Tbl = removevars(Tbl, varslist);
end


function Tbl = mySelect(Tbl, times, selbetween)
if selbetween
    times = timerange(times(1), times(2));
end
Tbl = Tbl(times,:);
Tbl.Properties.Events = Tbl.Properties.Events(times,:);
end


function evwinds = groupbytime(evtime, boundtime, timeBegin)
evtime = sort(evtime);
if isempty(evtime)
    evwinds = repmat(datetime(NaT,'TimeZone',timeBegin.TimeZone), 0,2);
else
    evbound = seconds(diff(evtime)) > boundtime;
    evend = evtime([evbound; false]);
    evstart = evtime([false; evbound]);
    evstart = [evtime(1); evstart];
    evend = [evend; evtime(end)];
    evwinds = [evstart, evend];
end
end


function [trng, tbls, chrej] = nancheck(trng, tbls, hAXs)
chrej = [];
tnan = []; cnan = {}; 
for SFi = 1:width(tbls)
    TBLi = tbls{SFi}; CHi = (TBLi.Properties.VariableNames);
    inan = sum(isnan(TBLi.Variables),2);
    tnan = [tnan; TBLi.Time(find(inan))];
    inan = sum(isnan(TBLi.Variables),1);
    cnan = [cnan, CHi(find(inan))];
end
cnan = unique(cnan);
cnans = cellfun(@(c) [c,'|'], cnan, 'UniformOutput', false);
cnans = [cnans{:}];
if ~isempty(tnan)
    for hAX = hAXs'
        plot(hAX, tnan, zeros(size(tnan)), '*r');
    end
    %xlim(trng)
    trimsel = questdlg(...
        ['Selection contains missing/nan values on channel(s) ',cnans], ...
        'Reselect due to nan values', ...
        'Start after nan', 'End before nan', 'Reject channels / Proceed', ...
        'Reject channels / Proceed');
    if isempty(trimsel)
        error('Selection must be made.')
    elseif strcmp(trimsel, 'Start after nan')
        trng(1) = max(tnan) + seconds(.001);
    elseif strcmp(trimsel, 'End before nan')
        trng(2) = min(tnan) - seconds(.001);
    end
    if strcmp(trimsel, 'Reject channels / Proceed')
        channelIndexReject = listdlg("PromptString","REJECT Channel(s)", ...
            "ListString",cnan, "SelectionMode","multiple"); 
        chrej = cnan(channelIndexReject); chrej = string(chrej);
    else
        tbls = cellfun(@(T) mySelect(T, trng, true), ...
            tbls, 'UniformOutput',false);
    end
end
end


function [trng, com, tblMainOut, tblsOut] = ...
    usrconfirm(trng, com, tblMainIn, tblsIn, hAXs, colr)

doAgain = true;
while doAgain
[trng, com] = VariableCountIntervalSelector(trng', com);
trng = trng'; 
[ucom,~,uInd] = unique(com);
tblMainOut = cell(length(ucom),1); tblsOut = cell(length(ucom),length(tblsIn));
hPat = [];
for ti = 1:width(trng)
    trng_ = trng(:,ti);
    hPat = [hPat; plottrng(trng_, colr, hAXs)]; % indicate on plot
    tblMainOut{uInd(ti)} = [tblMainOut{uInd(ti)}; mySelect(tblMainIn, trng_, true)];
    for SFi = 1:width(tblsOut)
        tblsOut{uInd(ti),SFi} = [tblsOut{uInd(ti),SFi}; mySelect(tblsIn{SFi}, trng_, true)];
    end
end

pause(.01); drawnow; pause(.01);
doAgainQst = questdlg('Accept?', 'Accept Time Range', ...
    'Yes', 'Try Again', 'Stop Program', 'Yes');
if strcmp(doAgainQst, 'Yes')
    doAgain = false;
elseif strcmp(doAgainQst, 'Try Again')
    delete(hPat);
else
    error('Terminated by user.')
end
end

for ic = 1:height(tblMainOut)
    %tblMain{ic,1} = myRetime(tblMain{ic,1}, SamplingFreq, nan);
    %tblMain{ic,1} = retime(tblMain{ic,1}, 'regular', 'nearest', 'SampleRate',SamplingFreq);
    tblMainOut{ic,1}.Properties.Description = ucom(ic);
    for SFi = 1:width(tblsOut)
        tblsOut{ic,SFi}.Properties.Description = ucom(ic);
    end
end
tblsOut = [tblMainOut, tblsOut];
if isequal(size(tblMainOut), [1,1])
    tblMainOut = tblMainOut{1};
end

end


function hPat = plottrng(trng, colr, hAXs)
hPat = [];
xrng = trng;
xrng = [xrng, xrng]'; xrng = xrng(:)';
xrng = [xrng, xrng(1)]; 
for hAX = hAXs'
    yrng = ylim(hAX); 
    yrng = [yrng, fliplr(yrng)];
    yrng = [yrng, yrng(1)];
    hPat = [hPat; patch(hAX, xrng, yrng, colr, 'FaceAlpha',.5)];
end
end


function plt = plotstem(vals, names)
plt = stem(vals); 
xticks(1:length(vals)); xticklabels(names);
end


function [img, txt] = plotgrid(vals, names, nrow, ncol)
%{
for n = 1:length(names)
    names{n} = [num2str(n),': ',names{n}];
end
%}
valgrid = nan(nrow, ncol);
for idx = 1:min(nrow*ncol, length(vals))
    valgrid(idx) = vals(idx);
end
img = imagesc(valgrid); colorbar;
hold on;
[X,Y] = meshgrid(1:ncol, 1:nrow);
names = names(1:min((ncol*nrow), length(names)));
X = X(:); X = X(1:length(names));
Y = Y(:); Y = Y(1:length(names));
txt = text(X,Y, names, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontWeight', 'bold', ...
    'Color',[.8 0 0]);
end


function saveasmultiple(fig, filename)
saveas(fig, filename, 'fig'); % original matlab figure
saveas(fig, filename, 'png'); % preview
end