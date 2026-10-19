%% Synthetic inputs for offline_PhaseDetect
% Generates two 40-second, 1 kHz sine-wave test signals and runs the
% offline phase detector on each one with plots and all outputs enabled.

clearvars;
close all;
clc;

%% Sampling and time vectors
SamplingFreq = 1000;                 % Hz
duration = 40;                       % seconds
numSamples = SamplingFreq * duration;
t = (0:numSamples-1) / SamplingFreq; % seconds; 0 through 39.999 s

%% Piecewise-constant-frequency sine wave: 20, 22, 16, then 20 Hz
segmentDuration = 10;                % seconds per segment
segmentSamples = segmentDuration * SamplingFreq;
piecewiseFrequency = [20, 22, 16, 20];
signalPiecewise = zeros(1, numSamples);
phaseAtSegmentStart = 0;

for k = 1:numel(piecewiseFrequency)
    idx = (k-1)*segmentSamples + (1:segmentSamples);
    localTime = (0:segmentSamples-1) / SamplingFreq;
    f = piecewiseFrequency(k);
    signalPiecewise(idx) = sin(phaseAtSegmentStart + 2*pi*f*localTime);
    phaseAtSegmentStart = phaseAtSegmentStart + 2*pi*f*segmentDuration;
end

%% Triangular frequency sweep, completing five 16-to-22-to-16 Hz cycles
% Each sweep cycle consists of a linear ramp from 16 to 22 Hz and back to
% 16 Hz. Integrating the instantaneous frequency gives continuous phase.
lowFrequency = 16;                   % Hz
highFrequency = 22;                  % Hz
numSweepCycles = 10;
sweepPhase = mod(t / duration * numSweepCycles, 1);
sweepFrequency = lowFrequency + (highFrequency-lowFrequency) .* ...
    (1 - abs(2*sweepPhase - 1));
signalSweep = sin(2*pi*cumtrapz(t, sweepFrequency));

%% Optional overview plot of the generated signals
figure('Name', 'Synthetic phase-detection inputs');
tiledlayout(2,1);
nexttile;
plot(t, signalPiecewise);
grid on;
xlabel('Time (s)'); ylabel('Amplitude');
title('Piecewise sine: 20, 22, 16, 20 Hz (10 s each)');
nexttile;
plot(t, signalSweep);
grid on;
xlabel('Time (s)'); ylabel('Amplitude');
title('Sine with five triangular frequency sweeps between 16 and 22 Hz');

%% Run detector on both inputs
% artDur = -1 disables artifact replacement. stepsize = 0 disables
% adaptive AR coefficient updates. showplots = true requests the full
% diagnostic plots from offline_PhaseDetect; dodebug = false avoids the
% per-sample debug animation.
StimTrainRec = false(size(t));
channelName = 'Synthetic';
PhaseOfInterest = 0;
FreqRange = [13, 30];
ARwin = 1000;
ARlen = 40;
predWin = 800;
artDur = -1;
packetLength = 10;
stepsize = 10;
donorm = 10;
doresample = false;
useIIR = true;
showplots = true;
dodebug = false;
PDSusefuture = true;

[piecewise.phAll, piecewise.phEst, piecewise.frAll, piecewise.frEst, ...
    piecewise.TtoStim, piecewise.phStim, piecewise.W, ...
    piecewise.dur, piecewise.numCycle, piecewise.numMissing, ...
    piecewise.numExtra] = offline_PhaseDetect( ...
    signalPiecewise, SamplingFreq, StimTrainRec, t, channelName, ...
    PhaseOfInterest, FreqRange, ARwin, ARlen, predWin, artDur, ...
    packetLength, stepsize, donorm, doresample, useIIR, showplots, ...
    dodebug, PDSusefuture);
%%
[sweep.phAll, sweep.phEst, sweep.frAll, sweep.frEst, ...
    sweep.TtoStim, sweep.phStim, sweep.W, sweep.dur, ...
    sweep.numCycle, sweep.numMissing, sweep.numExtra] = ...
    offline_PhaseDetect(signalSweep, SamplingFreq, StimTrainRec, t, ...
    channelName, PhaseOfInterest, FreqRange, ARwin, ARlen, predWin, ...
    artDur, packetLength, stepsize, donorm, doresample, useIIR, ...
    showplots, dodebug, PDSusefuture);

disp('Synthetic phase-detection runs complete.');
