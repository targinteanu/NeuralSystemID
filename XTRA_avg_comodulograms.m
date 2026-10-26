root = '/Users/torenarginteanu/Documents/Anderson Lab/XTRA/Detailed Heat Map Plots/2026-10 surrogate/PD25N009';
subjects = {'LS52', 'LS53', 'LS54'};
groups = {'Baseline no Notch', 'NBack no Notch'};
allData = cell(2, 3);

% Extract CData from each figure
for g = 1:2
    for s = 1:3
        figPath = fullfile(root, ...
            [subjects{s} ' Figures'], ...
            groups{g}, ...
            ['PAC_zMI_' subjects{s} '.fig']);

        fig = openfig(figPath);%, 'invisible');
        imgs = findobj(fig, 'Type', 'image');

        if isempty(imgs)
            close(fig);
            error('No image object found in %s', figPath);
        end

        allData{g, s} = get(imgs(1), 'CData');
        f_ph = get(imgs(1), 'XData');
        f_Amp = get(imgs(1), 'YData');
        close(fig);
    end
end

%% Average each group, checking matrix sizes
avg = cell(1, 2);
for g = 1:2
    referenceSize = size(allData{g, 1});

    for s = 2:3
        if ~isequal(size(allData{g, s}), referenceSize)
            error('Matrix sizes differ within group %s', groups{g});
        end
    end

    avg{1,g} = (double(allData{g, 1}) + ...
              double(allData{g, 2}) + ...
              double(allData{g, 3})) / 3;
end

if ~isequal(size(avg{1,1}), size(avg{1,2}))
    error('Baseline and NBack matrix sizes differ');
end

%{
% Save the averages and source matrices
save(fullfile(root, 'PAC_zMI_averages.mat'), ...
    'avg', 'allData', 'subjects', 'groups');

writematrix(avg{1}, fullfile(root, 'PAC_zMI_baseline_average.csv'));
writematrix(avg{2}, fullfile(root, 'PAC_zMI_nback_average.csv'));
%}

% Use a shared color scale for the two heatmaps
colorLimits = [min([avg{1,1}(:); avg{1,2}(:)]), ...
               max([avg{1,1}(:); avg{1,2}(:)])];
%%
for g = 1:2
    fig = figure;%('Visible', 'off', 'Color', 'w');

    img = imagesc(f_ph, f_Amp, avg{1,g});
    img.Parent.YDir = 'normal';
    axis xy;
    axis tight;
    axis square;
    colorbar;
    caxis(colorLimits);

    title(sprintf('Mean PAC zMI: %s (n=3)', groups{g}), ...
          'Interpreter', 'none');
    xlabel('Phase (Hz)');
    ylabel('Amplitude (Hz)');

    outputName = sprintf('PAC_zMI_%s_average.png', ...
        lower(strrep(groups{g}, ' ', '_')));

    exportgraphics(fig, fullfile(root, outputName), 'Resolution', 180);
    %close(fig);
end