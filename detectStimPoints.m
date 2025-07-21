function stim_points = detectStimPoints(patient_data, baseline_data, threshold)
% detectStimPoints - Identify timepoints of stimulation artifacts
% original by Jack in testing4_15_JC.m / edited by Toren 
%
% Inputs:
%   patient_data   - neural data matrix (channels x timepoints)
%   baseline_data  - baseline segment of data (channels x timepoints)
%   threshold      - fraction of channels that must exceed threshold (e.g., 0.5)
%
% Output:
%   stim_points - vector of time indices where artifact is detected

    % === Derivative across time ===
    d_data = diff(patient_data, 1, 2);            % channels x (timepoints - 1)
    d_data_abs = abs(d_data);

    % === Derivative of baseline ===
    d_baseline = diff(baseline_data, 1, 2);
    d_baseline_abs = abs(d_baseline);

    % === Baseline stats per channel ===
    baseline_mean = mean(d_baseline_abs, 2);      % [channels x 1]
    baseline_std  = std(d_baseline_abs, 0, 2);    % [channels x 1]

    stim_mask = d_data_abs > (baseline_mean + 10 * baseline_std);  % logical mask

    % === Proportion of channels above threshold per timepoint ===
    stim_proportion = mean(stim_mask, 1);         % 1 x (timepoints - 1)

    % === Timepoints where > threshold fraction of channels fire ===
    stim_points = find(stim_proportion > threshold) - 1;  % shift for diff()
end