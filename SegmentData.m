%% user selects folder; necessary files are pulled 
disp('Please select subject folder...')
folder = uigetdir; 
[~,pName] = fileparts(folder)
NSfiles = dir([folder,filesep,'*.ns*']);
NEVfiles = dir([folder,filesep,'*.nev']);
ElecXL = [dir([folder,filesep,'electrode_data.xls*']); ...
          dir([folder,filesep,pName,filesep,'electrode_data.xls*'])];

%% get electrode data 
disp('Reading electrode localization data...')
if isempty(ElecXL)
    warning('No electrode localization data found!')
end

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

clear x y z XYZ xyz fTbl

%% get blackrock data into a time/event table 
disp('Accessing blackrock data...')

timeBegin = datetime([inf inf inf], 'TimeZone','UTC'); 
timeEnd = datetime([-inf -inf -inf], 'TimeZone','UTC');

% NS (continuous data) files 
disp('  - continuous data...')
tbls = {}; SampleRates = [];
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
                tbls{SFi} = [tbls{SFi}; fTbl];
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
disp('  - event data...')
NEVtbl = [];
for f = NEVfiles
    try
    fEV = openNEV([f.folder,filesep,f.name], 'nosave');
    EVtbl = nev2table(fEV);
    timeBegin = min(timeBegin, min(EVtbl.Time));
    timeEnd = max(timeEnd, max(EVtbl.Time));
    NEVtbl = [NEVtbl; EVtbl];
    catch ME
        warning(ME.message);
    end
end

chnames = {};
for SFi = 1:width(tbls)
    % pair NS with NEV
    tbls{SFi}.Properties.Events = [tbls{SFi}.Properties.Events; NEVtbl];
    % get all channel names 
    chnames = [chnames, tbls{SFi}.Properties.VariableNames];
end
chnames = unique(chnames); chnames = string(chnames);

clear fNS fEV fTbl EVtbl c r cname

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

%% detect triggers 
disp('Detecting triggers...')

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
            tbls{SFi} = NStbl;
        end
    end
end

% group triggers by time 
trigboundtime = 10; % sec (defines a burst of triggers)
trigtime = trigtbl.Time; % all trigger types!
if isempty(trigtime)
    trigwinds = repmat(datetime(NaT,'TimeZone',timeBegin.TimeZone), 0,2);
else
trigbound = seconds(diff(trigtime)) > trigboundtime;
trigend = trigtime([trigbound; false]); 
trigstart = trigtime([false; trigbound]);
trigstart = [trigtime(1); trigstart];
trigend = [trigend; trigtime(end)];
trigwinds = [trigstart, trigend];
end

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

clear NStbl t trig lbl trigname trigtbl_ trignames_ trigname_

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

%% determine the baseline 
% can any of this be changed to use MainTable? 
disp('Detecting baseline...')

% assess file start/end, packet loss, etc
discont = [];
for SFi = 1:width(tbls)
    evtbl = tbls{SFi}.Properties.Events;
    for disctype = ["Start of file", "End of file", "SerialDigitalIO"]
        discevt = contains(evtbl.EventLabels, disctype);
        discont = [discont; evtbl.Time(discevt)];
    end
    disctype = "Data Loss";
        discevt = contains(evtbl.EventLabels, disctype);
        discont = [discont; evtbl.EventEnds(discevt)];
    inan = isnan(tbls{SFi}.Variables); inan = sum(inan,2);
    discont = [discont; tbls{SFi}.Time(find(inan))];
end

clear disctype discevt inan evtbl

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

clear trigstart_ trigend_ t1 t2 disci d

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

clear t1 t2 numOL io trng

%% user confirms baseline 
disp('Please confirm baseline condition...')
[~,iwind] = min(OLwinds);
fig1 = figure('Units','normalized', 'Position',[.05,.05,.9,.9]); 
[~,hAXs] = myStackedPlot(MainTable, ...
    [channelNameRec, channelNameStim, channelNameInspect]);
trngBaseline = inputdlg({'Start', 'End'}, 'BASELINE Segment', ...
    1, string(candwinds(iwind,:)));
trngBaseline = datetime(trngBaseline, 'TimeZone',candwinds.TimeZone);
tblBaselineMain = mySelect(MainTable, trngBaseline, true);
tblsBaseline = cellfun(@(T) mySelect(T, trngBaseline, true), ...
    tbls, 'UniformOutput',false);
tblsBaseline = [{tblBaselineMain}, tblsBaseline];

% nan check
[trngBaseline,tblsBaseline] = nancheck(trngBaseline,tblsBaseline,hAXs);
tblBaselineMain = tblsBaseline{1};

% indicate on plot
plottrng(trngBaseline, [0,0,1], hAXs);

pause(.01); drawnow; pause(.01);

%% show channel stats
% altered version of tblChannelSummary.m
disp('Calculating baseline characteristics...')
BaselineData = tblBaselineMain.Variables;
varnames = tblBaselineMain.Properties.VariableNames; 
vardescs = tblBaselineMain.Properties.VariableDescriptions;
vardescs = cellfun(@(d) [d(1:3),'...',d((end-20):end)], ...
    vardescs, 'UniformOutput',false);

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

clear showas

%% inspect and reject noisy channels 
disp('Please inspect and specify noisy channels to reject...')

% final inspection 
channelIndexInspect2 = listdlg("PromptString","Inspect Channel(s)", ...
    "ListString",chnames, "SelectionMode","multiple"); 
channelNameInspect2 = chnames(channelIndexInspect2);
fig4 = figure('Units','normalized', 'Position',[.05,.05,.9,.9]); 
myStackedPlot(MainTable, channelNameInspect2);

% rejection 
channelIndexReject = listdlg("PromptString","REJECT Channel(s)", ...
    "ListString",chnames, "SelectionMode","multiple"); 
channelNameReject = chnames(channelIndexReject);
MainTable = removevars(MainTable, channelNameReject);
tblBaselineMain = removevars(tblBaselineMain, channelNameReject);
tbls = cellfun(@(tbl) removevars(tbl, channelNameReject), ...
    tbls, 'UniformOutput',false);
tblsBaseline = cellfun(@(tbl) removevars(tbl, channelNameReject), ...
    tblsBaseline, 'UniformOutput',false);

%% user confirms marked stimulation 
disp('Please confirm stimulation triggers...')
figure(fig1);
[trngTrig, comTrig] = VariableCountIntervalSelector(trigwinds, trignames);
trngTrig = trngTrig'; 
[comTrig,~,uInd] = unique(comTrig);
tblTrigMain = cell(length(comTrig),1); tblsTrig = cell(length(comTrig),length(tbls));
for itrig = 1:width(trngTrig)
    trngTrig_ = trngTrig(:,itrig);
    plottrng(trngTrig_, [1,0,0], hAXs); % indicate on plot
    tblTrigMain{uInd(itrig)} = [tblTrigMain{uInd(itrig)}; mySelect(MainTable, trngTrig_, true)];
    for SFi = 1:width(tblsTrig)
        tblsTrig{uInd(itrig),SFi} = [tblsTrig{uInd(itrig),SFi}; mySelect(tbls{SFi}, trngTrig_, true)];
    end
end
for ic = 1:height(tblTrigMain)
    %tblTrigMain{ic,1} = myRetime(tblTrigMain{ic,1}, SamplingFreq, nan);
    %tblTrigMain{ic,1} = retime(tblTrigMain{ic,1}, 'regular', 'nearest', 'SampleRate',SamplingFreq);
    tblTrigMain{ic,1}.Properties.Description = comTrig(ic);
    for SFi = 1:width(tblsTrig)
        tblsTrig{ic,SFi}.Properties.Description = comTrig(ic);
    end
end
tblsTrig = [tblTrigMain, tblsTrig];
if isequal(size(tblTrigMain), [1,1])
    tblTrigMain = tblTrigMain{1};
end

pause(.01); drawnow; pause(.01);

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
if isempty(stimtime)
    stimwinds = repmat(datetime(NaT,'TimeZone',timeBegin.TimeZone), 0,2);
else
stimbound = seconds(diff(stimtime)) > stimboundtime;
stimend = stimtime([stimbound; false]); 
stimstart = stimtime([false; stimbound]);
stimstart = [stimtime(1); stimstart];
stimend = [stimend; stimtime(end)];
stimwinds = [stimstart, stimend];
end

% user confirms 
figure(fig1);
[trngStimNoTrig, comStim] = VariableCountIntervalSelector(stimwinds, ...
    repmat("Presumed Stim", height(stimwinds), 1));
trngStimNoTrig = trngStimNoTrig';
[comStim,~,uInd] = unique(comStim);
tblStimNoTrigMain = cell(length(comStim),1); 
tblsStimNoTrig = cell(length(comStim),length(tbls));
for istim = 1:width(trngStimNoTrig)
    trngStim_ = trngStimNoTrig(:,istim);
    plottrng(trngStim_, [1,0,1], hAXs); % indicate on plot
    tblStimNoTrigMain{uInd(istim)} = ...
        [tblStimNoTrigMain{uInd(istim)}; mySelect(MainTable, trngStim_, true)];
    for SFi = 1:width(tblsStimNoTrig)
        tblsStimNoTrig{uInd(istim),SFi} = ...
            [tblsStimNoTrig{uInd(istim),SFi}; ...
            mySelect(tbls{SFi}, trngStim_, true)];
    end
end
for ic = 1:height(tblStimNoTrigMain)
    %tblStimNoTrigMain{ic,1} = myRetime(tblStimNoTrigMain{ic,1}, SamplingFreq, nan);
    %tblStimNoTrigMain{ic,1} = retime(tblStimNoTrigMain{ic,1}, 'regular', 'nearest', 'SampleRate',SamplingFreq);
    tblStimNoTrigMain{ic,1}.Properties.Description = comStim(ic);
    for SFi = 1:width(tblsStimNoTrig)
        tblsStimNoTrig{ic,SFi}.Properties.Description = comStim(ic);
    end
end
tblsStimNoTrig = [tblStimNoTrigMain, tblsStimNoTrig];
if isequal(size(tblStimNoTrigMain), [1,1])
    tblStimNoTrigMain = tblStimNoTrigMain{1};
end

pause(.01); drawnow; pause(.01);

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
if isempty(srltime)
    srlwinds = repmat(datetime(NaT,'TimeZone',timeBegin.TimeZone), 0,2);
else
srlbound = seconds(diff(srltime)) > srlboundtime;
srlend = srltime([srlbound; false]); 
srlstart = srltime([false; srlbound]);
srlstart = [srltime(1); srlstart];
srlend = [srlend; srltime(end)];
srlwinds = [srlstart, srlend];
end
% marked stim 
if isempty(srltimeMarkedStim)
    srlwindsMarkedStim = repmat(datetime(NaT,'TimeZone',timeBegin.TimeZone), 0,2);
else
srlbound = seconds(diff(srltimeMarkedStim)) > srlboundtime;
srlend = srltimeMarkedStim([srlbound; false]); 
srlstart = srltimeMarkedStim([false; srlbound]);
srlstart = [srltimeMarkedStim(1); srlstart];
srlend = [srlend; srltimeMarkedStim(end)];
srlwindsMarkedStim = [srlstart, srlend];
end
% unmarked stim 
if isempty(srltimeUnmarkedStim)
    srlwindsUnmarkedStim = repmat(datetime(NaT,'TimeZone',timeBegin.TimeZone), 0,2);
else
srlbound = seconds(diff(srltimeUnmarkedStim)) > srlboundtime;
srlend = srltimeUnmarkedStim([srlbound; false]); 
srlstart = srltimeUnmarkedStim([false; srlbound]);
srlstart = [srltimeUnmarkedStim(1); srlstart];
srlend = [srlend; srltimeUnmarkedStim(end)];
srlwindsUnmarkedStim = [srlstart, srlend];
end
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
[trngSrl, comSrl] = VariableCountIntervalSelector(srlwinds, srlnames);
trngSrl = trngSrl';
[comSrl,~,uInd] = unique(comSrl);
tblSrlMain = cell(length(comSrl),1); 
tblsSrl = cell(length(comSrl),length(tbls));
for isrl = 1:width(trngSrl)
    trngSrl_ = trngSrl(:,isrl);
    plottrng(trngSrl_, [0,1,0], hAXs); % indicate on plot
    tblSrlMain{uInd(isrl)} = ...
        [tblSrlMain{uInd(isrl)}; mySelect(MainTable, trngSrl_, true)];
    for SFi = 1:width(tblsSrl)
        tblsSrl{uInd(isrl),SFi} = ...
            [tblsSrl{uInd(isrl),SFi}; ...
            mySelect(tbls{SFi}, trngSrl_, true)];
    end
end
for ic = 1:height(tblSrlMain)
    %tblSrlMain{ic,1} = myRetime(tblSrlMain{ic,1}, SamplingFreq, nan);
    %tblSrlMain{ic,1} = retime(tblSrlMain{ic,1}, 'regular', 'nearest', 'SampleRate',SamplingFreq);
    tblSrlMain{ic,1}.Properties.Description = comSrl(ic);
    for SFi = 1:width(tblsSrl)
        tblsSrl{ic,SFi}.Properties.Description = comSrl(ic);
    end
end
tblsSrl = [tblSrlMain, tblsSrl];
if isequal(size(tblSrlMain), [1,1])
    tblSrlMain = tblSrlMain{1};
end

pause(.01); drawnow; pause(.01);

clear srlname srltbl_ srlnames_ srlname_ trngSrl_

%% segment misc 
disp('Please label any other data to save...')
% user confirms 
figure(fig1);
[trngMisc, comMisc] = VariableCountIntervalSelector(...
    repmat(datetime(NaT,'TimeZone',timeBegin.TimeZone), 1,2), ...
    "Segment any miscellaneous time periods, or close this window to cancel.");
trngMisc = trngMisc';
[comMisc,~,uInd] = unique(comMisc);
tblMiscMain = cell(length(comMisc),1); 
tblsMisc = cell(length(comMisc),length(tbls));
for imisc = 1:width(trngMisc)
    trngMisc_ = trngMisc(:,imisc);
    plottrng(trngMisc_, [1,1,0], hAXs); % indicate on plot
    tblMiscMain{uInd(imisc)} = ...
        [tblMiscMain{uInd(imisc)}; mySelect(MainTable, trngMisc_, true)];
    for SFi = 1:width(tblsMisc)
        tblsMisc{uInd(imisc),SFi} = ...
            [tblsMisc{uInd(imisc),SFi}; ...
            mySelect(tbls{SFi}, trngMisc_, true)];
    end
end
for ic = 1:height(tblMiscMain)
    %tblMiscMain{ic,1} = myRetime(tblMiscMain{ic,1}, SamplingFreq, nan);
    %tblMiscMain{ic,1} = retime(tblMiscMain{ic,1}, 'regular', 'nearest', 'SampleRate',SamplingFreq);
    tblMiscMain{ic,1}.Properties.Description = comMisc(ic);
    for SFi = 1:width(tblsMisc)
        tblsMisc{ic,SFi}.Properties.Description = comMisc(ic);
    end
end
tblsMisc = [tblMiscMain, tblsMisc];
if isequal(size(tblMiscMain), [1,1])
    tblMiscMain = tblMiscMain{1};
end

pause(.01); drawnow; pause(.01);

clear trngMisc_

%% save all
disp('Saving segmented data...')

thisfilename = mfilename("fullpath");
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
    "tblsMisc", "tblsSrl", "tblsStimNoTrig", "tblsTrig", "tblsBaseline", ...
    "channelNameRec", "channelNameStim", "channelNameTrig", "channelNameReject", ...
    "trngMisc", "trngSrl", "trngStimNoTrig", "trngTrig", "trngBaseline", ...
    "-v7.3");

% save fig1 and fig4 
disp('  - figure 1:')
saveasmultiple(fig1, [svname,'_TimePlot1']);
disp('  - figure 4:')
saveasmultiple(fig4, [svname,'_TimePlot2']);
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

function Tbl = mySelect(Tbl, times, selbetween)
if selbetween
    times = timerange(times(1), times(2));
end
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
        tbls = cellfun(@(T) mySelect(T, trngBasleline, true), ...
            tbls, 'UniformOutput',false);
    end
end
end

function plottrng(trng, colr, hAXs)
xrng = trng;
xrng = [xrng, xrng]'; xrng = xrng(:)';
xrng = [xrng, xrng(1)]; 
for hAX = hAXs'
    yrng = ylim(hAX); 
    yrng = [yrng, fliplr(yrng)];
    yrng = [yrng, yrng(1)];
    patch(hAX, xrng, yrng, colr, 'FaceAlpha',.5);
end
end

function plt = plotstem(vals, names)
plt = stem(vals); 
xticks(1:length(vals)); xticklabels(names);
end

function [img, txt] = plotgrid(vals, names, nrow, ncol)
for n = 1:length(names)
    names{n} = [num2str(n),': ',names{n}];
end
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