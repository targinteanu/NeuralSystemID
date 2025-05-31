%% user selects folder; necessary files are pulled 

% e-phys data
folder = uigetdir('', 'Select Electrophysiological Data Folder.'); 
[~,pName] = fileparts(folder)
NS2files = dir([folder,filesep,'*.ns2']);
NEVfiles = dir([folder,filesep,'*.nev']);

% imaging data
folderImg = uigetdir(folder, 'Select Imaging Data Folder.'); 
elec_acpc_fr_file = dir([folderImg,filesep,'*elec_acpc_fr.mat']);
elec_fsavg_frs_file = dir([folderImg,filesep,'*elec_fsavg_frs.mat']);
pial_lh_file = [folderImg,filesep,'freesurfer',filesep,'surf',filesep,'lh.pial.T1'];
pial_rh_file = [folderImg,filesep,'freesurfer',filesep,'surf',filesep,'rh.pial.T1'];
for flist = ["elec_acpc_fr_file", "elec_fsavg_frs_file"]
    if isempty(eval(flist))
        error(flist+" list is empty!");
    end
    if length(eval(flist)) > 1
        fnames = {flist.name};
        sel = listdlg("PromptString","Select "+flist, "SelectionMode","single", "ListString",fnames);
        eval(flist+" = "+flist+"(sel);");
    end
    eval(flist+" = fullfile("+flist+".folder,"+flist+".name);");
end
clear sel flist fnames

%% get blackrock data into a time/event table 
NS2tbl = [];
for f = NS2files
    fNS = openNSx([f.folder,filesep,f.name], 'uV');
    fTbl = ns2timetable(fNS);
    NS2tbl = [NS2tbl; fTbl];
end

NEVtime = [];
for f = NEVfiles
    try
        fEV = openNEV([f.folder,filesep,f.name]);
        ft0 = datetime(fEV.MetaTags.DateTime);
        ftRel = fEV.Data.SerialDigitalIO.TimeStampSec;
        fTime = ft0 + seconds(ftRel);
        NEVtime = [NEVtime; fTime];
    catch ME
        warning(ME.message)
    end
end
if ~isempty(NEVtime)
    NS2tbl.Properties.Events = eventtable(NEVtime, ...
        EventLabels="Serial Digital IO");
end

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

%% Beta Power -> color map
BetaColor = BetaPower(1:(end-1))'; % exclude ainp1 
BetaColor = normalize(BetaColor, "range"); 
BetaColor = [1,1,0].*BetaColor; % yellow

%% main plotting - edited 

ft_defaults

pial_lh = ft_read_headshape(pial_lh_file);
pial_lh.coordsys = 'acpc';
pial_rh = ft_read_headshape(pial_rh_file);
pial_rh.coordsys = 'acpc';

load(elec_acpc_fr_file);
figure; 
ft_plot_mesh(pial_lh);
ft_plot_mesh(pial_rh);
ax = gca;
ax.Children(1).FaceColor = [.8 .8 .8];
ax.Children(2).FaceColor = [.8 .8 .8];
el=ft_plot_sens(elec_acpc_fr);
el.ZData = el.ZData + 5;
el.CData = BetaColor;
el.SizeData = 50;
el.Marker = "o";
el.MarkerFaceColor = 'flat';
el.MarkerEdgeColor = 'r';
% view([66.1890, 39.71]);
material shiny;
lighting gouraud;
camlight;
camlight headlight;

%% main plotting - Jack
%{
ft_defaults

load_ACPC_FR_mesh(elec_acpc_fr_file, pial_lh_file, BetaPower);

load_FSAVG_FRS_mesh(elec_fsavg_frs_file, pial_lh_file, BetaPower);

min_val = min(BetaPower);
max_val = max(BetaPower);
BetaPowerNormalized = (BetaPower - min_val) / (max_val - min_val);
BetaPowerNormalized(BetaPowerNormalized >= 1) = 0.999;

load_ACPC_FR_mesh(elec_acpc_fr_file, pial_lh_file, BetaPowerNormalized);
%}

%% helper function(s) 

function load_ACPC_FR_mesh(acpc_path, pial_lh_path, data)

    % === Step 1: Load Electrode Positions and Mesh ===
    elec_acpc_fr = load(acpc_path).elec_acpc_fr;
    pial_lh = ft_read_headshape(pial_lh_path);
    pial_lh.coordsys = 'acpc';

    pial_lh = ft_convert_units(pial_lh, 'mm');
    elec_acpc_fr = ft_convert_units(elec_acpc_fr, 'mm');

    % === Step 2: Identify reference channels ===
    ref_idx = isnan(data);          % Logical array: true if it's reference
    data_clean = data;              % Copy data
    data_clean(ref_idx) = 0;         % Set NaNs to 0 for interpolation

    % === Step 3: Build source structure ===
    source = struct();
    source.pos = elec_acpc_fr.chanpos;
    source.pow = data_clean;
    source.dimord = 'pos';
    source.coordsys = 'acpc';

    % === Step 4: Interpolate electrode data ===
    cfg = [];
    cfg.parameter = 'pow';
    cfg.interpmethod = 'sphere_weighteddistance';
    cfg.sphereradius = 3;
    interp_source = ft_sourceinterpolate(cfg, source, pial_lh);

    % === Step 5: Plot the Interpolated Data ===
    cfg = [];
    cfg.funparameter = 'pow';
    cfg.funcolorlim = [0 1];            % normal 0-1 scaling
    cfg.method = 'surface';
    cfg.interpmethod = 'nearest';
    cfg.sphereradius = 3;
    cfg.camlight = 'no';

    ft_sourceplot(cfg, interp_source, pial_lh);

    view([-55 10]);
    material dull;
    lighting gouraud;
    camlight;

    % === Step 6: Normal color map (pure parula) ===
    colormap(parula);   % NO special red added

    % === Step 7: Plot electrodes ===
    % Plot normal electrodes
    elec_plot_colors = data_clean;   % Use interpolated data (0–1)
    ft_plot_sens(elec_acpc_fr, 'elecsize', 35, 'facecolor', elec_plot_colors, 'edgecolor', 'black');

    % === Step 8: Overlay reference electrodes manually ===
    if any(ref_idx)
        ref_positions = elec_acpc_fr.chanpos(ref_idx, :);
        hold on;
        scatter3(ref_positions(:,1), ref_positions(:,2), ref_positions(:,3), ...
            100, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        hold off;
    end

end



function load_FSAVG_FRS_mesh(fsavg_path, pial_lh_path, data)

    % === Step 1: Load Electrode Positions and Mesh ===
    elec_fsavg = load(fsavg_path).elec_fsavg_frs;
    pial_lh = ft_read_headshape('lh.pial');
    pial_lh.coordsys = 'fsaverage';

    pial_lh = ft_convert_units(pial_lh, 'mm');
    elec_fsavg = ft_convert_units(elec_fsavg, 'mm');

    % === Step 2: Prepare Electrode Data ===
    ref_idx = isnan(data);           % Find reference channels
    elec_values = data;

    % Replace NaNs with a very large dummy value
    dummy_value = 2;  % well above PLV 1
    elec_values(ref_idx) = dummy_value;  % mark reference channels

    % === Step 3: Build source structure ===
    source = struct();
    source.pos = elec_fsavg.chanpos;
    source.inside = true(size(source.pos,1),1);
    source.pow = elec_values;
    source.dimord = 'pos';
    source.coordsys = 'fsaverage';

    % === Step 4: Interpolate electrode data ===
    cfg = [];
    cfg.parameter = 'pow';
    cfg.interpmethod = 'sphere_weighteddistance';
    cfg.sphereradius = 8;
    interp_source = ft_sourceinterpolate(cfg, source, pial_lh);

    % === Step 5: After interpolation, restore reference channel to a clean value ===
    % Any interpolated value >1.5 is from dummy (reference channels)
    ref_mask = interp_source.pow > 1.5;
    interp_source.pow(ref_mask) = 1.001;  % set a standard "red" value

    % === Step 6: Plot brain surface ===
    cfg = [];
    cfg.funparameter = 'pow';
    cfg.funcolorlim = [0 1];  % PLVs go between 0 and 1
    cfg.method = 'surface';
    cfg.interpmethod = 'nearest';
    cfg.sphereradius = 3;
    cfg.camlight = 'no';

    ft_sourceplot(cfg, interp_source, pial_lh);

    view([-55 10]);
    material dull;
    lighting gouraud;
    camlight;

    % === Step 7: Build custom colormap ===
    base_cmap = parula(256);
    red_color = [1 0 0];  % pure red
    cmap_with_red = [base_cmap; red_color];
    colormap(gca, cmap_with_red);

    % === Step 8: Also plot electrodes ===
    elec_plot_colors = data;
    elec_plot_colors(ref_idx) = 1.001;  % manually set reference channels to red
    ft_plot_sens(elec_fsavg, 'elecsize', 35, 'facecolor', elec_plot_colors, 'edgecolor', 'black');

end