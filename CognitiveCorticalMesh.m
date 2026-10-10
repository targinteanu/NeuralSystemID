%% obtain segmented, artifact-free data 

% file selection
[fn,fp] = uigetfile('*.mat', 'Choose Artifact-Free Data File');
load(fullfile(fp,fn));
[~,fn,fe] = fileparts(fn);
ArtRemoveDone = contains(fn, '_ArtifactRemoveOffline');
QCcleaned = contains(fn,'QCcleaned');
if ArtRemoveDone
    fnOrig = split(fn, '_ArtifactRemoveOffline'); fnOrig = fnOrig{1};
    fpOrig = fp;
    while ~exist(fullfile(fpOrig,[fnOrig,fe]), 'file')
        fpOrig = fileparts(fpOrig); % try one folder out
    end
    load(fullfile(fpOrig,[fnOrig,fe]));
end

%% data processing

% clean baseline selection 
PDdata = tblsBaseline{1};
srate = PDdata.Properties.SampleRate;
if isnan(srate)
    srate = 1/median(seconds(diff(PDdata.Time)));
end

%% channel selection
channames = num2str((1:63)');
channames(channames==' ') = '0';
channames = "LS"+string(channames);
channames = channames';
for ch = channames
    if ~sum(strcmp(PDdata.Properties.VariableNames, ch))
        % fill missing channel with nan
        PDdata.(ch) = nan(height(PDdata),1);
    end
end
PDdata = PDdata(:,channames); % reorder

% common avg ECoG reref 
PDdata{:,:} = PDdata{:,:} - mean(PDdata{:,1:63}, 2, 'omitnan');

% interp nan values 
PDdata_nanmask = double(isnan(PDdata.Variables));
for ch = 1:width(PDdata)
    x = PDdata{:,ch};
    if sum(~isnan(x))
        if sum(isnan(x))
            iAll = (1:length(x))'; iValid = find(~PDdata_nanmask(:,ch));
            xValid = x(iValid);
            x = interp1(iValid,xValid,iAll,"linear","extrap");
            PDdata{:,ch} = x;
        end
    end
end

% signal processing 
%thetaPowerCortex = bandpower(PDdata.Variables, srate, [4,9]);
%thetaPowerCortex(isoutlier(sqrt(thetaPowerCortex))) = nan;
filtwts = fir1(1024, [4, 9]./(srate/2));
%PDdata = FilterTimetable(@(b,x) filtfilt(b,1,x), filtwts, PDdata);
for ch = 1:width(PDdata)
    disp(['Filtering channel ',num2str(ch),' of ',num2str(width(PDdata))])
    x = PDdata{:,ch};
    if sum(~isnan(x))
        PDdata{:,ch} = filtfilt(filtwts,1,x);
    end
end
PDdata = PDdata(:,1:63); PDdata_nanmask = PDdata_nanmask(:,1:63);
%PD_Phase_Data = instPhaseFreqTbl(PDdata);
PD_Channel_Names = PDdata.Properties.VariableNames;

ElecTbl = table('RowNames',PD_Channel_Names(1:63)); % cortical only

%% windowed theta power 
% Calculate windowed theta power
windowSize = 10000; % samples 
thetaPowerWindowed = PDdata.Variables.^2;
for ch = 1:width(PDdata)
    disp(['Calculating power: channel ',num2str(ch),' of ',num2str(width(PDdata))])
    x = thetaPowerWindowed(:,ch);
    if sum(~isnan(x))
        thetaPowerWindowed(:,ch) = movmean(envelope(x.^2), windowSize);
    end
end
PDdata_nanmask = movmean(PDdata_nanmask, windowSize);

%thetaPowerWindowed = movmean(envelope(PDdata.Variables).^2, windowSize);
thetaPowerWindowed = thetaPowerWindowed(round(windowSize/2):windowSize:end, :);
PDdata_nanmask = PDdata_nanmask(round(windowSize/2):windowSize:end, :);
PDdata_nanmask = PDdata_nanmask>0;
thetaPowerWindowed = 10*log10((thetaPowerWindowed)); % decibel (dB) scale 
if ~QCcleaned
    thetaPowerWindowed(isoutlier(thetaPowerWindowed)) = nan;
    % if QC already done, this step is not necessary?
end
thetaPowerWindowed(PDdata_nanmask) = nan;
thetaPowerWindowed = median(thetaPowerWindowed, 'omitnan');

ElecTbl.ThetaPowerWindowed = thetaPowerWindowed';
G = zeros(21,3);
G(:) = thetaPowerWindowed; 
R = mean(G,2, 'omitnan');

OLthresh = thetaPowerWindowed > median(thetaPowerWindowed, 'omitnan');
OLthresh = thetaPowerWindowed(OLthresh); 
io = isoutlier(OLthresh, 'mean'); OLthresh = min(OLthresh(io));
if ~QCcleaned
    % if QC already done, this step is not necessary?
    disp(['Removing ',num2str(sum(io)),' outlier channels.'])
    if isempty(OLthresh) || isnan(OLthresh)
        OLthresh = inf;
    end
    thetaPowerWindowed(thetaPowerWindowed >= OLthresh) = nan;
end

%%
%{
figure; 
G = zeros(21,3);
subplot(1,2,1); 
G(:) = log10(thetaPowerCortex); 
imagesc(G); colorbar;
title('bandpower')
subplot(1,2,2);
G(:) = log10(thetaPowerWindowed);
imagesc(G); colorbar;
title('median windowed envelope')
%}
thetaPowerOrig = thetaPowerWindowed;

%% obtain imaging data
ft_defaults
[ftver, ftpath] = ft_version;

% specify laterality 
lat = questdlg('Specify Electrode Laterality:', 'Laterality', ...
    'left','right','Cancel','Cancel');
if ~(strcmp(lat,'left')|strcmp(lat,'right'))
    error('Laterality must be specified.')
else
    if strcmp(lat,'left')
        lat = 1;
    else
        lat = 2;
    end
end

% specify file(s)
[fnElec, fpElec] = uigetfile('*_elec_acpc_f*.mat');
subjID = split(fnElec, '_elec_acpc'); subjID = subjID{1};
fnElecMNI = [subjID,'_elec_mni_frv.mat'];

% look for imaging files 
fpPialT1 = fpElec;
fnPialT1 = fullfile(fpPialT1,'lh.pial.T1');
if ~exist(fnPialT1, 'file')
    fpPialT1 = fullfile(fpElec,'freesurfer','surf');
    fnPialT1 = fullfile(fpPialT1,'lh.pial.T1');
end
meshACPCl = ft_read_headshape(fnPialT1);
fnPialT1 = fullfile(fpPialT1,'rh.pial.T1');
if ~exist(fnPialT1, 'file')
    fpPialT1 = fullfile(fpElec,'freesurfer','surf');
    fnPialT1 = fullfile(fpPialT1,'rh.pial.T1');
end
meshACPCr = ft_read_headshape(fnPialT1);

% load files
elecACPC = load(fullfile(fpElec, fnElec)).elec_acpc_fr;
elecMNI = load(fullfile(fpElec, fnElecMNI)).elec_mni_frv;
meshMNIl = load([ftpath,filesep,'template/anatomy/surface_pial_left.mat']).mesh;
meshMNIr = load([ftpath,filesep,'template/anatomy/surface_pial_right.mat']).mesh;
meshACPC = {meshACPCl, meshACPCr};
meshMNI = {meshMNIl, meshMNIr};

%% patient specific acpc coords

ElecTbl.acpcX = elecACPC.chanpos(:,1); 
ElecTbl.acpcY = elecACPC.chanpos(:,2); 
ElecTbl.acpcZ = elecACPC.chanpos(:,3); 

for ch = 1:length(thetaPowerWindowed)
    if isnan(thetaPowerWindowed(ch))
        dr = elecACPC.chanpos(ch,:) - elecACPC.chanpos;
        dl = rms(dr,2); L = 1./dl';
        L = L([1:(ch-1), (ch+1):end]); theta = thetaPowerOrig([1:(ch-1), (ch+1):end]);
        thetaPowerWindowed(ch) = sum(theta.*L, 'omitnan')/sum(L, 'omitnan');
    end
end

figACPC = figure('Units','normalized', 'Position',[.05,.05,.65,.9], 'Theme','light'); 
ThetaPowerInterp = showbrainmesh(elecACPC, meshACPC{lat}, ...
    thetaPowerWindowed, 'acpc', meshACPC{3-lat});%, isnan(thetaPowerOrig));

pos = mean(elecACPC.chanpos);
[az,el] = cart2sph(pos(1),pos(2),pos(3));

camlight;
%{
if lat(1)=='r'
    view([150,5]);
else
    view([0,90]);
end
%}
view(az*180/pi + 90, el*180/pi);
sgtitle('Theta Power, individual ACPC', 'FontSize',20)
camlight headlight;

%% freesurfer mni coords 

ElecTbl.mniX = elecMNI.chanpos(:,1); 
ElecTbl.mniY = elecMNI.chanpos(:,2); 
ElecTbl.mniZ = elecMNI.chanpos(:,3); 

figMNI = figure('Units','normalized', 'Position',[.3,.05,.65,.9], 'Theme','light'); 
ThetaPowerInterp = showbrainmesh(elecMNI, meshMNI{lat}, ...
    thetaPowerWindowed, 'mni', meshMNI{3-lat});%, isnan(thetaPowerOrig));

pos = mean(elecMNI.chanpos);
[az,el] = cart2sph(pos(1),pos(2),pos(3));

camlight;
%{
if lat(1)=='r'
    view([150,5]);
else
    view([0,90]);
end
%}
view(az*180/pi + 90, el*180/pi);
sgtitle('Theta Power, MNI space / template brain', 'FontSize',20)
camlight headlight;

%% saving 
doSave = questdlg('Save?');
if strcmp(doSave, 'Yes')
writetable(ElecTbl, fullfile(fpElec,'BrainHeatmapCoords.xlsx'));
writematrix(R, fullfile(fpElec,'R_row_vector.xlsx'));
saveas(figACPC, fullfile(fpElec,'BrainHeatmapACPC'), 'fig');
saveas(figMNI, fullfile(fpElec,'BrainHeatmapMNI'), 'fig');
saveas(figMNI, fullfile(fpElec,'BrainHeatmapMNI'), 'png');
saveas(figACPC, fullfile(fpElec,'BrainHeatmapACPC'), 'png');
end

%% helper functions 

function [interp_source, srcpltcfg] = ...
    showbrainmesh(elecpos, brainmesh, data, coordsys, brainmesh2, markmissing)
    sphereradius = 4;
    if nargin < 6
        markmissing = false(size(elecpos.chanpos,1),1);
        if nargin < 5
            brainmesh2 = [];
        end
    end

    % === Step 1: Load Electrode Positions and Mesh ===
    brainmesh.coordsys = coordsys;

    brainmesh = ft_convert_units(brainmesh, 'mm');
    elecpos = ft_convert_units(elecpos, 'mm');

    %{
    % === Step 2: Identify reference channels ===
    ref_idx = isnan(data);          % Logical array: true if it's reference
    data_clean = data;              % Copy data
    data_clean(ref_idx) = 0;         % Set NaNs to 0 for interpolation
    %}

    % === Step 3: Build source structure ===
    source = struct();
    source.pos = elecpos.chanpos;
    source.pow = data;
    source.dimord = 'pos';
    source.coordsys = coordsys;
    source.unit = 'mm';

    % === Step 4: Interpolate electrode data ===
    cfg = [];
    cfg.parameter = 'pow';
    cfg.interpmethod = 'sphere_weighteddistance';
    cfg.sphereradius = sphereradius;
    interp_source = ft_sourceinterpolate(cfg, source, brainmesh);

    % === Step 5: Plot the Interpolated Data ===
    cfg = [];
    cfg.funparameter = 'pow';
    %cfg.funcolorlim = [0 1];            % normal 0-1 scaling
    cfg.funcolorlim = [min(data), max(data)]; % normalize
    cfg.method = 'surface';
    cfg.interpmethod = 'nearest';
    cfg.sphereradius = sphereradius;
    cfg.camlight = 'no';

    srcpltcfg = ft_sourceplot(cfg, interp_source, brainmesh);

    %view([-55 10]);
    ax = gca;
    ax.Children(2).FaceColor = [.8 .8 .8]; % brain mesh
    ax.Children(2).FaceAlpha = .6; % brain mesh
    material shiny;
    lighting gouraud;
    %camlight;
    %view([150, 5]);

    fig = gcf; 
    fig.Children(2).FontSize = 16; % colorbar 
    fig.Children(2).Label.String = 'Power (dB)';
    fig.Children(2).Label.FontSize = 18;

    % === Step 6: Normal color map (pure parula) ===
    colormap(parula);   % NO special red added

    % === Step 7: Plot electrodes ===
    % Plot normal electrodes
    elec_plot_colors = data;   % Use interpolated data (0–1)
    %ft_plot_sens(elec_acpc_fr, 'elecsize', 9, 'facecolor', elec_plot_colors, 'edgecolor', 'black');
    ft_plot_sens(elecpos, 'elecsize', 15, 'facecolor', 'red', 'edgecolor', 'black');

    % indicate electrodes with missing/omitted data 
    hold on;
    plot3(...
        elecpos.chanpos(markmissing,1), ...
        elecpos.chanpos(markmissing,2), ...
        elecpos.chanpos(markmissing,3)+.1, ...
        'xk', 'MarkerSize',15, 'LineWidth',3);

    if ~isempty(brainmesh2)
        hold on; ft_plot_mesh(brainmesh2);
    end
    ax.Children(1).FaceColor = [.8 .8 .8]; % brain mesh 2
    ax.Children(1).FaceAlpha = .6; % brain mesh 2

    %{
    % === Step 8: Overlay reference electrodes manually ===
    if any(ref_idx)
        ref_positions = elec_acpc_fr.chanpos(ref_idx, :);
        hold on;
        scatter3(ref_positions(:,1), ref_positions(:,2), ref_positions(:,3), ...
            100, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        hold off;
    end
    %}

end
