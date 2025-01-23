%% user selects folder; necessary files are pulled 
folder = uigetdir; 
[~,pName] = fileparts(folder)
NS2files = dir([folder,filesep,'*.ns2']);
NEVfiles = dir([folder,filesep,'*.nev']);

%% get blackrock data into a time/event table 
NS2tbl = [];
for f = NS2files
    fNS = openNSx([f.folder,filesep,f.name], 'uV');
    fTbl = ns2timetable(fNS);
    NS2tbl = [NS2tbl; fTbl];
end

NEVtime = [];
for f = NEVfiles
    fEV = openNEV([f.folder,filesep,f.name]);
    ft0 = datetime(fEV.MetaTags.DateTime);
    ftRel = fEV.Data.SerialDigitalIO.TimeStampSec;
    fTime = ft0 + seconds(ftRel);
    NEVtime = [NEVtime; fTime];
end
NS2tbl.Properties.Events = eventtable(NEVtime, ...
    EventLabels="Serial Digital IO");

clear f fNS fEV fTbl ft0 ftRel fTime NEVtime

%% get details on stimulus and recording 

tblDur = @(T) T.Time(end) - T.Time(1);

% assume stim pulse train is ainp1 - can change this to user selection 
channelIndexStimTrain = find(contains(NS2tbl.Properties.VariableNames, 'AINP1'));
if isempty(channelIndexStimTrain)
    StimTrig = false(height(NS2tbl), 1);
else
    channelIndexStimTrain = channelIndexStimTrain(1);
    StimTrig = diff(NS2tbl{:,channelIndexStimTrain}) > 1e3;
    StimTrig = [false; StimTrig];
end
StimTrigTime = NS2tbl.Time(StimTrig);
%NS2tbl = [NS2tbl, table(StimTrig)];

% data before first stim 
Stim1 = find(StimTrig);
if ~isempty(Stim1)
    StimEnd = Stim1(end); StimEnd = min(height(NS2tbl), StimEnd+5000);
    Stim1 = Stim1(1);
    PreStimEnd = max(1, Stim1-5000);
else
    PreStimEnd = height(NS2tbl);
    StimEnd = 1;
end
PreStimEndTime = NS2tbl.Time(PreStimEnd);
StimEndTime = NS2tbl.Time(StimEnd);
tblPreStim = NS2tbl(1:PreStimEnd,:);
disp(['Pre Stim: ',...
    char(tblDur(tblPreStim)),' (',num2str(PreStimEnd),' samples) out of '...
    char(tblDur(NS2tbl)),' (',num2str(height(NS2tbl)),' samples)'])

% summary channel data
[BetaPower, SD, ~, ~, fig1] = tblChannelSummary(tblPreStim, [13, 30]);
subplot(4,1,1); ylabel('Beta Band Power'); title('Channels Summary Data - Pre Stim');
sgtitle(pName);
isOut = isoutlier(NS2tbl, 'mean');

%% inspect or reject channels 

% user selects rec and stim channels 
channelIndexRec = listdlg("PromptString","Recording Channel(s)", ...
    "ListString",NS2tbl.Properties.VariableNames, "SelectionMode","multiple"); 
channelIndexStim = listdlg("PromptString","Stimulus Channel(s)", ...
    "ListString",NS2tbl.Properties.VariableNames, "SelectionMode","multiple"); 
channelIndexInspect = listdlg("PromptString","Inspect Channel(s)", ...
    "ListString",NS2tbl.Properties.VariableNames, "SelectionMode","multiple"); 
stimChanName = '';
for ch = channelIndexStim
    CH = NS2tbl.Properties.VariableNames{channelIndexStim};
    CH = CH(CH > 0); % get rid of white space
    stimChanName = [stimChanName,' and ',CH];
end
stimChanName = stimChanName(6:end);

chanSelInds = unique([channelIndexRec, channelIndexStim, ...
    channelIndexInspect, channelIndexStimTrain]);
fig2 = figure; myStackedPlot(NS2tbl(:,chanSelInds), [], isOut(:,chanSelInds)); 
sgtitle(pName);

% user selects bad channels 
channelIndexRem = listdlg("PromptString","Remove Channel(s)", ...
    "ListString",NS2tbl.Properties.VariableNames, "SelectionMode","multiple"); 
channelIndexKeep = true(1,width(NS2tbl)); 
channelIndexKeep(channelIndexRem) = false;
if ~isempty(channelIndexStimTrain)
    channelIndexKeep(channelIndexStimTrain) = false;
end
chanSelInds = chanSelInds(channelIndexKeep(chanSelInds));
tblSel = NS2tbl(:,channelIndexKeep); 
tblSel.Properties.Events = [tblSel.Properties.Events; ...
    eventtable(StimTrigTime, ...
    EventLabels="Stimulus on channel(s) "+string(stimChanName))];

%% select time range 
selTime = @(T, t1, t2) T( ((T.Time >= t1) & (T.Time <= t2)) , :);

tSel = inputdlg(...
    {'Start Time:', 'End Time:'}, ...
    'Select Time Range', ...
    1, ...
    {char(tblSel.Time(1)), char(tblSel.Time(end))});

if ~isempty(tSel)
    tSel = cellfun(@datetime, tSel);
    tblSel = selTime(tblSel, tSel(1), tSel(2));
end

%% assign data 

% quiet baseline #1 
BL1EndTime = PreStimEndTime;
useEvt = true;
try
    evt = tblSel.Properties.Events;
    evtAct = strcmp(evt.EventLabels, "Serial Digital IO");
    evtAct = evt(evtAct, :);
    evtTime = evtAct.Time;
    if ~isempty(evtTime)
        BL1EndTime = min(min(evtTime)-seconds(5), BL1EndTime);
    else
        useEvt = false;
    end
catch ME
    warning(['Events were not considered due to issue: ',ME.message])
    useEvt = false;
end
tblBL1 = selTime(tblSel, tblSel.Time(1), BL1EndTime);
tblSel = selTime(tblSel, BL1EndTime, tblSel.Time(end));

% activity with no stim
if useEvt
    evtTime = sort(evtTime);
    evtDur = median(diff(evtTime));
    tAct = [min(evtTime)-evtDur, max(evtTime)+evtDur];
    tAct(2) = min(tAct(2), PreStimEndTime);
    tblAct = selTime(tblSel, tAct(1), tAct(2));
    tblSel = selTime(tblSel, tAct(2), tblSel.Time(end));
else
    tblAct = tblSel([],:); % empty 
end

% recorded stimulus 
if StimEndTime >= PreStimEndTime
    tblStim = selTime(tblSel, PreStimEndTime, StimEndTime);
    tblSel = selTime(tblSel, StimEndTime, tblSel.Time(end));
else
    tblStim = tblSel([],:); % empty
end

% try to find instances of DBS 
% (stim not recorded)
% simultaneous artifact in most channels = DBS 
isOutSel = isoutlier(tblSel, 'mean');
figure; plot(sum(isOutSel,2)/size(isOutSel,2)); grid on;
dbsTrig = sum(isOutSel,2)/size(isOutSel,2) > .1;
dbs1 = find(dbsTrig);
if ~isempty(dbs1)
    dbsEnd = dbs1(end); dbsEnd = min(height(tblStim), dbsEnd+5000);
    dbs1 = dbs1(1);
    dbs1 = max(1, dbs1-5000);
    tblDBS = tblSel(dbs1:dbsEnd,:);

    % quiet baseline #2 
    if dbs1 > 1
        tblBL2 = tblSel(1:dbs1,:); 
    else
        tblBL2 = tblSel([],:); % empty
    end
else
    tblBL2 = tblSel;
    tblDBS = tblSel([],:); % empty
end

DataTimeTables = {...
    'Baseline 1',            trimEvents(tblBL1),  tblBL1.Time(1),  tblBL1.Time(end); ...
    'Activity Without Stim', trimEvents(tblAct),  tblAct.Time(1),  tblAct.Time(end); ...
    'Recorded Stim',         trimEvents(tblStim), tblStim.Time(1), tblStim.Time(end); ...
    'Baseline 2',            trimEvents(tblBL2),  tblBL2.Time(1),  tblBL2.Time(end); ...
    'DBS',                   trimEvents(tblDBS),  tblDBS.Time(1),  tblDBS.Time(end)}
for h = 1:height(DataTimeTables)
    figure('Units','normalized', 'Position',...
        [(h-1)/height(DataTimeTables),0,1/height(DataTimeTables),1]);
    myStackedPlot(DataTimeTables{h,2}(:,chanSelInds), [], []); % chanSelInds needs to be changed!
    sgtitle({pName; DataTimeTables{h,1}});
end

%% saving 
saveas(fig1, [folder,filesep,pName,'_ChannelSummary'],'fig'); 
saveas(fig1, [folder,filesep,pName,'_ChannelSummary'],'png'); 
saveas(fig2, [folder,filesep,pName,'_Data'],'fig'); 
saveas(fig2, [folder,filesep,pName,'_Data'],'png'); 
save([folder,filesep,pName,'_DataTimeTables.mat'], 'DataTimeTables');

%% helper
function tbl = trimEvents(tbl)
evt = tbl.Properties.Events;
evtTime = evt.Time;
t1 = tbl.Time(1); t2 = tbl.Time(end);
trimind = (evtTime >= t1) & (evtTime <= t2);
evt = evt(trimind, :);
tbl.Properties.Events = evt;
end