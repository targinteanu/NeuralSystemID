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
NS2tbl.Properties.Events = NEVtime;

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
NS2tbl = [NS2tbl, table(StimTrig)];

% data before first stim 
Stim1 = find(StimTrig);
if ~isempty(Stim1)
    Stim1 = Stim1(1);
    PreStimEnd = max(0, Stim1-5000);
else
    PreStimEnd = height(NS2tbl);
end
tblPreStim = NS2tbl(1:PreStimEnd,:);
disp(['Pre Stim: ',...
    char(tblDur(tblPreStim)),' (',num2str(PreStimEnd),' samples) out of '...
    char(tblDur(NS2tbl)),' (',num2str(height(NS2tbl)),' samples)'])

% summary channel data
% might need to change this if multiple NS2 files with break in between,
% i.e. non-uniform sampling 
BetaPower = bandpower(tblPreStim.Variables, tblPreStim.Properties.SampleRate, [13, 30]);
figure; 
subplot(2,1,1); stem(BetaPower); 
ylabel('Beta Band Power'); title('Channels Summary Data (Pre-Stim)');
xticks(1:width(BetaPower)); xticklabels(NS2tbl.Properties.VariableNames);
subplot(2,1,2); boxplot(tblPreStim.Variables, 'PlotStyle','compact', 'Symbol','.');
ylabel('Box Plot'); xlabel('Channel Name');
xticks(1:width(BetaPower)); xticklabels(NS2tbl.Properties.VariableNames);
sgtitle(pName);

%% inspect or reject channels 

% user selects rec and stim channels 
channelIndexRec = listdlg("PromptString","Recording Channel(s)", ...
    "ListString",NS2tbl.Properties.VariableNames, "SelectionMode","multiple"); 
channelIndexStim = listdlg("PromptString","Stimulus Channel(s)", ...
    "ListString",NS2tbl.Properties.VariableNames, "SelectionMode","multiple"); 
channelIndexInspect = listdlg("PromptString","Inspect Channel(s)", ...
    "ListString",NS2tbl.Properties.VariableNames, "SelectionMode","multiple"); 

chanSelInds = unique([channelIndexRec, channelIndexStim, ...
    channelIndexInspect, channelIndexStimTrain]);
figure; myStackedPlot(NS2tbl(:,chanSelInds)); grid on; sgtitle(pName);

% user selects bad channels 
channelIndexRem = listdlg("PromptString","Remove Channel(s)", ...
    "ListString",NS2tbl.Properties.VariableNames, "SelectionMode","multiple"); 
channelIndexKeep = true(1,width(NS2tbl)); 
channelIndexKeep(channelIndexRem) = false;
tblSel = NS2tbl(:,channelIndexKeep); 