%% obtain segmented, artifact-free data 
%%%%%% 60 hz notch 
% PAC settings
phaseFreqs = 3:0.25:9;
ampFreqs = 30:2:200;
phaseBandWidthHz = 1;
ampBandWidthHz = 10;
pacUseFullRecording = false;
pacWinSec = 10; % make this 10s
pacNBin = 18;

% file selection
[fn,fp] = uigetfile('*SegmentData*.mat', 'Choose Artifact-Free Data File');
load(fullfile(fp,fn));
[~,fn,fe] = fileparts(fn);
ArtRemoveDone = contains(fn, '_ArtifactRemoveOffline');
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
miscIdx = min(2, size(tblsMisc, 1));
PDdata = tblsMisc{miscIdx, 1};
srate = PDdata.Properties.SampleRate;
if isnan(srate)
    srate = 1/median(seconds(diff(PDdata.Time)));
end
if max(ampFreqs) >= srate/2
    error(['Sampling rate (',num2str(srate), ' Hz) is too low for ', ...
    'requested amplitude band up to ', num2str(max(ampFreqs)), ' Hz.']);
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
PDraw = PDdata{:,1:63}; % preserve broadband rereferenced data for PAC

%{
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
PDdata = PDdata(:,1:63);
%PD_Phase_Data = instPhaseFreqTbl(PDdata);
%}
PD_Channel_Names = PDdata.Properties.VariableNames;

%% choose channels for PAC comodulograms
[selIdx, ok] = listdlg( ...
    'PromptString', 'Select channel(s) for PAC comodulograms:', ...
    'SelectionMode', 'multiple', ...
    'ListString', cellstr(PD_Channel_Names), ...
    'ListSize', [300, 400]);
if ~ok || isempty(selIdx)
    error('No channel selected for PAC comodulograms.')
end
pacChannelNames = PD_Channel_Names(selIdx);

%% phase-amplitude coupling (PAC) comodulograms
disp('Computing PAC comodulograms: phase 3-9 Hz, amplitude 30-200 Hz...')
PACResults = struct();
PACFigures = cell(numel(selIdx), 1);
for ii = 1:numel(selIdx)
    ch = selIdx(ii);
    chName = pacChannelNames{ii};
    fprintf('Computing PAC comodulogram for %s (%d of %d)\n', chName, ii, numel(selIdx));
    pacMat = compute_pac_comodulogram( ...
        PDraw(:,ch), srate, phaseFreqs, ampFreqs, ...
        phaseBandWidthHz, ampBandWidthHz, pacUseFullRecording, pacWinSec, pacNBin);

    PACResults(ii).channelIndex = ch;
    PACResults(ii).channelName = chName;
    PACResults(ii).phaseFreqs = phaseFreqs;
    PACResults(ii).ampFreqs = ampFreqs;
    PACResults(ii).pacMat = pacMat;

    figPAC = figure('Units','normalized', 'Position',[.1,.1,.6,.7], ...
        'Name', ['PAC_', chName], 'NumberTitle', 'off');
    imagesc(phaseFreqs, ampFreqs, pacMat.');
    axis xy;
    xlabel('Frequency for phase (Hz)');
    ylabel('Frequency for amplitude (Hz)');
    title(['Phase-Amplitude Coupling: ', chName], 'Interpreter','none');
    colormap(turbo);
    cb = colorbar;
    cb.Label.String = 'Tort MI';
    cb.Label.FontSize = 12;
    set(gca, 'FontSize', 12, 'LineWidth', 1.2);
    PACFigures{ii} = figPAC;
end

ElecTbl = table('RowNames',PD_Channel_Names(1:63)); % cortical only

%% windowed theta power 
% Calculate windowed theta power
windowSize = 1000; % samples 
thetaPowerWindowed = PDdata.Variables.^2;
for ch = 1:width(PDdata)
    disp(['Calculating power: channel ',num2str(ch),' of ',num2str(width(PDdata))])
    x = thetaPowerWindowed(:,ch);
    if sum(~isnan(x))
        thetaPowerWindowed(:,ch) = movmean(envelope(x.^2), windowSize);
    end
end
%thetaPowerWindowed = movmean(envelope(PDdata.Variables).^2, windowSize);
thetaPowerWindowed = thetaPowerWindowed(round(windowSize/2):windowSize:end, :);
thetaPowerWindowed = 10*log10((thetaPowerWindowed)); % decibel (dB) scale 
thetaPowerWindowed(isoutlier(thetaPowerWindowed)) = nan;
thetaPowerWindowed = median(thetaPowerWindowed, 'omitnan');

ElecTbl.ThetaPowerWindowed = thetaPowerWindowed';
G = zeros(21,3);
G(:) = thetaPowerWindowed; 
R = mean(G,2, 'omitnan');

OLthresh = thetaPowerWindowed > median(thetaPowerWindowed, 'omitnan');
OLthresh = thetaPowerWindowed(OLthresh); 
io = isoutlier(OLthresh, 'mean'); OLthresh = min(OLthresh(io));
disp(['Removing ',num2str(sum(io)),' outlier channels.'])
if isempty(OLthresh) || isnan(OLthresh)
    OLthresh = inf;
end
thetaPowerWindowed(thetaPowerWindowed >= OLthresh) = nan;

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
end

% specify file(s)
[fnElec, fpElec] = uigetfile('*_elec_acpc_f*.mat');
subjID = split(fnElec, '_elec_acpc'); subjID = subjID{1};
fnElecMNI = [subjID,'_elec_mni_frv.mat'];

% look for imaging files 
fpPialT1 = fpElec;
fnPialT1 = fullfile(fpPialT1,[lat(1),'h.pial.T1']);
if ~exist(fnPialT1, 'file')
    fpPialT1 = fullfile(fpElec,'freesurfer','surf');
    fnPialT1 = fullfile(fpPialT1,[lat(1),'h.pial.T1']);
end

% load files
elecACPC = load(fullfile(fpElec, fnElec)).elec_acpc_fr;
elecMNI = load(fullfile(fpElec, fnElecMNI)).elec_mni_frv;
meshACPC = ft_read_headshape(fnPialT1);
meshMNI = load([ftpath,filesep,'template/anatomy/surface_pial_',lat,'.mat']).mesh;

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

figACPC = figure('Units','normalized', 'Position',[.05,.05,.65,.9]); 
ThetaPowerInterp = load_ACPC_FR_mesh(elecACPC, meshACPC, ...
    thetaPowerWindowed, 'acpc', isnan(thetaPowerOrig));
title('Theta Power, individual ACPC', 'FontSize',20)

camlight;
if lat(1)=='r'
    view([150,5]);
else
    view([0,90]);
end
camlight headlight;

%% freesurfer mni coords 

ElecTbl.mniX = elecMNI.chanpos(:,1); 
ElecTbl.mniY = elecMNI.chanpos(:,2); 
ElecTbl.mniZ = elecMNI.chanpos(:,3); 

figMNI = figure('Units','normalized', 'Position',[.3,.05,.65,.9]); 
ThetaPowerInterp = load_ACPC_FR_mesh(elecMNI, meshMNI, ...
    thetaPowerWindowed, 'mni', isnan(thetaPowerOrig));
title('Theta Power, MNI space / template brain', 'FontSize',20)

camlight;
if lat(1)=='r'
    view([150,5]);
else
    view([0,90]);
end
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
save(fullfile(fpElec, 'PACComodulograms.mat'), 'PACResults', 'phaseFreqs', 'ampFreqs', 'phaseBandWidthHz', 'ampBandWidthHz', 'pacUseFullRecording', 'pacWinSec', 'pacNBin');
for ii = 1:numel(PACResults)
    saveas(PACFigures{ii}, fullfile(fpElec, ['PAC_', PACResults(ii).channelName]), 'fig');
    saveas(PACFigures{ii}, fullfile(fpElec, ['PAC_', PACResults(ii).channelName]), 'png');
end
end

%% helper functions 

function pacMat = compute_pac_comodulogram(x, srate, phaseFreqs, ampFreqs, phaseBandWidthHz, ampBandWidthHz, useFullRecording, winSec, nbin)
% Compute a PAC comodulogram per channel.
% Default mode uses full-recording MI for each frequency pair.
    nPhase = numel(phaseFreqs);
    nAmp = numel(ampFreqs);
    pacMat = nan(nPhase, nAmp);

    for ip = 1:nPhase
        phaseBand = [max(phaseFreqs(ip) - phaseBandWidthHz/2, 0.1), phaseFreqs(ip) + phaseBandWidthHz/2];
        phaseBand(2) = min(phaseBand(2), srate/2 - 0.1);
        if phaseBand(2) <= phaseBand(1)
            continue;
        end

        for ia = 1:nAmp
            ampBand = [max(ampFreqs(ia) - ampBandWidthHz/2, 0.1), ampFreqs(ia) + ampBandWidthHz/2];
            ampBand(2) = min(ampBand(2), srate/2 - 0.1);
            if ampBand(2) <= ampBand(1)
                continue;
            end

            if useFullRecording
                pacMat(ip, ia) = compute_full_pac_tort(x, srate, phaseBand, ampBand, nbin);
            else
                pacMat(ip, ia) = compute_windowed_pac_tort(x, srate, phaseBand, ampBand, winSec, nbin);
            end
        end
    end
end

function mi = compute_full_pac_tort(x, srate, phaseBandHz, ampBandHz, nbin)
% Compute Tort MI using the full recording (single estimate).
    nSamp = size(x,1);
    xc = x(:);
    mi = nan;
    if all(isnan(xc))
        return;
    end

    valid = ~isnan(xc);
    if sum(valid) < 10
        return;
    end
    if any(~valid)
        ti = (1:nSamp)';
        xc(~valid) = interp1(ti(valid), xc(valid), ti(~valid), 'linear', 'extrap');
    end

    xc = apply_line_notches(xc, srate);
    bPhase = fir1(1024, phaseBandHz./(srate/2));
    bAmp = fir1(1024, ampBandHz./(srate/2));
    xPhase = filtfilt(bPhase, 1, xc);
    xAmp = filtfilt(bAmp, 1, xc);

    % calcPAC performs Hilbert transform and trims center 80% to reduce edge effects.
    mi = calcPAC(xPhase, xAmp, nbin, []);
end

function pacVals = compute_windowed_pac_tort(x, srate, phaseBandHz, ampBandHz, winSec, nbin)
% Compute median Tort MI across fixed-length windows for each channel.
    nCh = size(x,2);
    nSamp = size(x,1);
    winSamp = max(round(winSec*srate), 1);

    bPhase = fir1(1024, phaseBandHz./(srate/2));
    bAmp = fir1(1024, ampBandHz./(srate/2));

    pacVals = nan(1, nCh);
    for ch = 1:nCh
        xc = x(:,ch);
        if all(isnan(xc))
            continue;
        end

        valid = ~isnan(xc);
        if sum(valid) < 10
            continue;
        end
        if any(~valid)
            ti = (1:nSamp)';
            xc(~valid) = interp1(ti(valid), xc(valid), ti(~valid), 'linear', 'extrap');
        end

        xc = apply_line_notches(xc, srate);
        xPhase = filtfilt(bPhase, 1, xc);
        xAmp = filtfilt(bAmp, 1, xc);

        hPhi = angle(hilbert(xPhase));
        hAmp = abs(hilbert(xAmp));

        idxStart = 1:winSamp:nSamp;
        winMI = nan(1, numel(idxStart));
        for w = 1:numel(idxStart)
            i1 = idxStart(w);
            i2 = min(i1 + winSamp - 1, nSamp);
            if (i2 - i1 + 1) < max(round(0.5*winSamp), 10)
                continue;
            end
            winMI(w) = calcPAChelper(hPhi(i1:i2), hAmp(i1:i2), nbin, []);
        end

        pacVals(ch) = median(winMI, 'omitnan');
    end
end

function x = apply_line_notches(x, srate)
% Remove 60, 120, and 180 Hz line-noise harmonics.
    notchFreqs = [60, 120, 180];
    qFactor = 35;
    nyquist = srate/2;

    for notchFreq = notchFreqs
        if notchFreq < nyquist
            omega = 2*pi*notchFreq/srate;
            radius = 1 - pi*notchFreq/(qFactor*srate);
            bNotch = [1, -2*cos(omega), 1];
            aNotch = [1, -2*radius*cos(omega), radius^2];
            x = filtfilt(bNotch, aNotch, x);
        end
    end
end

function [interp_source, srcpltcfg] = ...
    load_ACPC_FR_mesh(elec_acpc_fr, pial_lh, data, coordsys, markmissing)
    sphereradius = 3;
    if nargin < 5
        markmissing = false(size(elec_acpc_fr.chanpos,1),1);
    end

    % === Step 1: Load Electrode Positions and Mesh ===
    %elec_acpc_fr = load(acpc_path).elec_acpc_fr;
    %pial_lh = ft_read_headshape(pial_lh_path);
    pial_lh.coordsys = coordsys;

    pial_lh = ft_convert_units(pial_lh, 'mm');
    elec_acpc_fr = ft_convert_units(elec_acpc_fr, 'mm');

    %{
    % === Step 2: Identify reference channels ===
    ref_idx = isnan(data);          % Logical array: true if it's reference
    data_clean = data;              % Copy data
    data_clean(ref_idx) = 0;         % Set NaNs to 0 for interpolation
    %}
    data_clean = data;

    % === Step 3: Build source structure ===
    source = struct();
    source.pos = elec_acpc_fr.chanpos;
    source.pow = data_clean;
    source.dimord = 'pos';
    source.coordsys = coordsys;
    source.unit = 'mm';

    % === Step 4: Interpolate electrode data ===
    cfg = [];
    cfg.parameter = 'pow';
    cfg.interpmethod = 'sphere_weighteddistance';
    cfg.sphereradius = sphereradius;
    interp_source = ft_sourceinterpolate(cfg, source, pial_lh);

    % === Step 5: Plot the Interpolated Data ===
    cfg = [];
    cfg.funparameter = 'pow';
    %cfg.funcolorlim = [0 1];            % normal 0-1 scaling
    cfg.funcolorlim = [min(data_clean), max(data_clean)]; % normalize
    cfg.method = 'surface';
    cfg.interpmethod = 'nearest';
    cfg.sphereradius = sphereradius;
    cfg.camlight = 'no';

    srcpltcfg = ft_sourceplot(cfg, interp_source, pial_lh);

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
    elec_plot_colors = data_clean;   % Use interpolated data (0–1)
    %ft_plot_sens(elec_acpc_fr, 'elecsize', 9, 'facecolor', elec_plot_colors, 'edgecolor', 'black');
    ft_plot_sens(elec_acpc_fr, 'elecsize', 15, 'facecolor', 'red', 'edgecolor', 'black');

    % indicate electrodes with missing/omitted data 
    hold on;
    plot3(...
        elec_acpc_fr.chanpos(markmissing,1), ...
        elec_acpc_fr.chanpos(markmissing,2), ...
        elec_acpc_fr.chanpos(markmissing,3)+.1, ...
        'xk', 'MarkerSize',15, 'LineWidth',3);

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
