%% offline_PhaseDetect 
% 
% Toren Arginteanu, Johns Hopkins School of Medicine 
% 
% This script loads recorded brain recording data and performs the phase
% detection algorithm as if it were receiving the data one sample at a time
% in real time. 
% Then, it performs the same phase detection using the full data as a
% ground truth, and compares the results of the simulated "real time"
% against the ground truth. 
% For artifact removal, the script will attempt to determine when the
% actual stimuli in the recorded data occurred, using the trigger channel
% if available and using outlier detection to detect stimulus artifacts.
% Note that for real time purposes, the triggers being scheduled by the
% real-time phase detection algorithm would be used to time the artifact
% removal instead. 
% 
% An overview of the phase detection algorithm is as follows: ...

function [phAll, phEst, frAll, frEst, toStim, phStim, W, R, dur] = ...
    offline_PhaseDetect(dataOneChannel, SamplingFreq, StimTrainRec, t, channelName, ...
    PhaseOfInterest, FreqRange, ARwin, ARlen, predWin, artDur, packetLength, ...
    stepsize, donorm, doresample, showplots)

if height(dataOneChannel) > 1
    error('Data should be one channel only and horizontal.')
end

% signal to use default values if any arguments are not passed in
if nargin < 16
    showplots = false;
    if nargin < 15
        doresample = false;
        if nargin < 14
            donorm = false;
            if nargin < 13
                stepsize = [];
                if nargin < 12
                    packetLength = [];
                    if nargin < 11
                        artDur = [];
                        if nargin < 10
                            predWin = [];
                            if nargin < 9
                                ARlen = [];
                                if nargin < 8
                                    ARwin = [];
                                    if nargin < 7
                                        FreqRange = [];
                                        if nargin < 6
                                            PhaseOfInterest = [];
                                            if nargin < 5
                                                channelName = '';
                                                if nargin < 4
                                                    t = [];
                                                    if nargin < 3
                                                        StimTrainRec = [];
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Constants: 


% Simulate phase-dependent stimulation at this phase: 
PhaseOfInterest_default = 0; % radians; i.e. 0 for peak, pi for trough stimulation

% frequency range: 
loco_default = 13; hico_default = 30; % low and high cutoff (Hz); e.g. 13-30 for beta band

% Artifact duration: 
artDur_default = 10; % artifact duration is extended by __ samples 

% AR model parameters: 
ARwin_default = 1000; % #samples of baseline to use to fit the AR model
ARlen_default = 10; % AR model order 

% AR-model prediction duration: 
predWin_default = 20; % #samples ahead to predict at each time step 

% length of packets of incomming data 
packetLength_default = 1; % samples; e.g. 1 for each sample entered individually 


% apply default values if necessary 
if isempty(FreqRange)
    hico = hico_default; loco = loco_default;
else
    hico = FreqRange(2); loco = FreqRange(1);
end
varnames = ["PhaseOfInterest", "artDur", "ARwin", "ARlen", "predWin", "packetLength"];
for v = varnames
    if isempty(eval(v))
        eval(v+" = "+v+"_default;");
    end
end
clear v varnames
if isempty(stepsize)
    if donorm
        stepsize = 1e-3;
    else
        % use data autocorrelation to get step size
        stepsize = LearnrateEst(dataOneChannel,ARlen);
        stepsize = stepsize/1000; % for additional safety
    end
end


samplerat = 1;
if doresample
    samplerat = 10; % downsample by this factor (must be integer). TO DO: replace with something based on nyquist rate!
end

%% Part A: Setup 
% Load and preprocess brain recording data from a file. 
% Setup the signals, filter, and AR model that will be used to simulate the
% real-time process. 

if nargin < 1
    [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
        channelName, channelIndex, channelIndexStim, channelNames]...
        = getRecordedData_NS();
else
    if isempty(StimTrainRec)
        StimTrainRec = false(size(dataOneChannel));
    end
    if isempty(t)
        t = 1:length(dataOneChannel); t = (t-1)/SamplingFreq;
    end
    if isempty(channelName)
        channelName = ' ';
    end
end

dataOneChannel = dataOneChannel - mean(dataOneChannel);
dataOneChannelWithArtifact = dataOneChannel; 

% Get indexes of stimulus: 
% find when artifacts are believed to occur 
isOut = isoutlier(dataOneChannel, 'mean');
isArt = isOut | StimTrainRec; 
if artDur > 0
    isArt = movsum(isArt, artDur) > 0;
end

%% A.2 Identify baseline and fit AR model
% In real time, this would be determined by the research team at some point
% during the procedure, preferably when the electrodes have been inserted
% into a stable location but the stimulus has not yet begun. 

% Find a clean baseline to train the AR model. 
artIndAll = isArt;
artIndAll(1) = true; artIndAll(end) = true;
artIndAll = find(artIndAll);
[~,baselineStartInd] = max(diff(artIndAll));
baselineEndInd = artIndAll(baselineStartInd+1); 
baselineStartInd = artIndAll(baselineStartInd); 
dataBaseline = dataOneChannel(baselineStartInd:baselineEndInd); 
baselineWin = (baselineEndInd-baselineStartInd) + samplerat*[-3,3]*ARwin; 
baselineWin = baselineWin/2; baselineWin = round(baselineWin); 
baselineWin(1) = max(1,baselineWin(1)); 
baselineWin(2) = min(length(dataBaseline),baselineWin(2));

% Train the AR model on unfiltered data. 
dataBaseline1 = dataBaseline(baselineWin(1):baselineWin(2));
ARmdl_unfilt = ar(iddata(dataBaseline1', [], 1/SamplingFreq), ARlen, 'yw');

%% A.3 Band-Pass Filtering setup 
% Get FIR filter weights, filter the signal, and train another AR model on
% filtered data. 

% filtering bound rules 
minfac         = 1;    % this many (lo)cutoff-freq cycles in filter
min_filtorder  = 15;   % minimum filter length

% filter order 
if loco>0
    filtord = minfac*fix(SamplingFreq/loco);
elseif hico>0
    filtord = minfac*fix(SamplingFreq/hico);
end
if filtord < min_filtorder
    filtord = min_filtorder;
end

% build filter and filtered baseline 
filtwts = fir1(filtord, [loco, hico]./(SamplingFreq/2));
dataBaseline = filtfilt(filtwts,1,dataBaseline);
%[dataBaseline, filtwts] = Myeegfilt(dataBaseline,SamplingFreq,loco,hico);
filtord = length(filtwts);
filtinit = zeros(filtord-1,1); % FIR filter Initial Condition
%filtdelay = ceil(filtord/2); % delay (#samples) caused by FIR filter
filtdelay = filtord;
dataBaseline1 = dataBaseline(baselineWin(1):baselineWin(2));

% downsample, if applicable
downsampledFreq = SamplingFreq; samplerat = 1;
if doresample
    dataBaseline1 = movmean(dataBaseline1, samplerat);
    dataBaseline1 = downsample(dataBaseline1, samplerat);
    downsampledFreq = SamplingFreq/samplerat;
end

% setup filtered AR model
disp(['Data size is ',num2str(length(dataBaseline1)/ARlen),' times variable size.'])
ARmdl_filt = ar(iddata(dataBaseline1', [], 1/downsampledFreq), ARlen, 'yw');
ARmdl_filt = ARmdl_filt.A;
ARmdl_filt_orig = ARmdl_filt; ARmdl_filt_new = ARmdl_filt; 

%% Part B: Real-Time Algorithm 
% Loop through each sample of data and perform the algorithm as if the data
% were being received in real time. 
toStim = false(size(dataOneChannel)); % intended stimulus trigger pulses 
phEst = nan(size(dataOneChannel)); % instantaneous phase estimate (rad)
frEst = nan(size(dataOneChannel)); % instantaneous frequency estimate (Hz)
i2nextStim_prev = Inf; % #samples to next stim pulse 
dataOneChannelFilt = zeros(size(dataOneChannel));
dur = nan(size(dataOneChannel));

% track AR mdl updates 
w = -ARmdl_filt(2:end)/ARmdl_filt(1);
w = fliplr(w);
W = nan(length(dataOneChannel), length(w)); 
w0 = w;

% transform model into eigenvalue domain for stability
r = roots([1, -fliplr(w)])';
rAmp = abs(r); rAng = angle(r);
stbl = max(rAmp) <= 1;
if stbl
    rAmpArg = 1./rAmp; rAmpArg = rAmpArg-1; 
    rAmpArg = -log(rAmpArg);
else
    warning('Initial AR model estimate is not stable.')
end
R = nan(length(dataOneChannel), length(r));
r0 = r;
stbl = false;

% calculate "actual" ground truth data that will be compared afterwards
%dataOneChannelFilt2 = Myeegfilt(dataOneChannel,SamplingFreq,loco,hico);
dataOneChannelFilt2 = filtfilt(filtwts,1,dataOneChannel);
[phAll, frAll] = instPhaseFreq(dataOneChannelFilt2, SamplingFreq);
frAll = min(frAll, hico); frAll = max(frAll, loco);

progTick = .05; prog = 0; % track progress

for tind = packetLength:packetLength:length(dataOneChannel)
    tic
    tinds = (tind-packetLength+1):tind;

    % track progress
    prog = prog + packetLength/length(dataOneChannel);
    if prog > progTick
        prog = prog - progTick; 
        disp(['Progress: ',num2str(100*tind/length(dataOneChannel)),'%']);
    end

    % Step 1: Remove Artifact 
    % Replace artifact-corrupted signal with AR model-forecasted data.  
    if artDur > 0
        for ind1 = tinds(isArt(tinds))
            ind0 = ind1 - ARlen;
            if ind0 > 0
                dataOneChannel(ind1) = myFastForecastAR(ARmdl_unfilt, ...
                    dataOneChannel(ind0:(ind1-1))', 1);
            end
        end
    end

    % Step 2: Filter 
    % Extract data within desired frequency band, e.g. beta
    [dataOneChannelFilt(tinds),filtinit] = filter(filtwts,1, ...
        dataOneChannel(tinds), filtinit);

    ind0 = tinds(1) - ARwin;
    if ind0 > 0
        dataPast = dataOneChannelFilt(ind0:tind)';
        dataPastOrig = dataPast;
        % downsample, if applicable 
        if doresample
            dataPast = flipud(dataPast);
            dataPast = dataPast(1:samplerat:end);
            dataPast = flipud(dataPast);
        end

        % Step 3: update AR coefficients 
        % Update weights using gradient descent
        if stepsize > 0
            r1 = r;
            w = -ARmdl_filt(2:end)/ARmdl_filt(1); 
            w = fliplr(w); w1 = w;
            x = dataPast((end-ARlen):(end-1));
            if (artDur <= 0) || ~isArt(tind)
                if stbl
                    [w, r, rAmpArg, rAng] = updateWtsStbl(...
                        w, r, rAmpArg, rAng, x, dataOneChannelFilt(tind));
                else
                    w = updateWts(w, x, dataOneChannelFilt(tind));
                    r = roots([1, -fliplr(w)]);
                end
                if max(abs(r)) < 1 % ensure stability
                    ARmdl_filt_new = [norm(w)/norm(w0), -fliplr(w)];
                else
                    w = w1; r = r1;
                end
            end
        end

        % Step 4: AR-based forecasting 
        % Predict some duration ahead using the AR model; this will be used
        % to pad the Hilbert transform
        dataFuture = myFastForecastAR(ARmdl_filt_new, dataPast, predWin);
        W(tind,:) = w; R(tind,:) = r;
        if stepsize > 0
            % prevent updated model from blowing up 
            if norm(dataFuture) <= 10*norm(dataPast) % set blowup threshold here
                ARmdl_filt = ARmdl_filt_new;
            else
                % revert 
                dataFuture = myFastForecastAR(ARmdl_filt, dataPast, predWin);
                W(tind,:) = w1; R(tind,:) = r1;
                if norm(dataFuture) > 100*norm(dataPast)
                    dataFuture = myFastForecastAR(ARmdl_filt_orig, dataPast, predWin);
                    W(tind,:) = w0; R(tind,:) = r1;
                    if norm(dataFuture) > norm(dataPast)
                        %keyboard
                    end
                end
            end
        end
        % upsample, if applicable 
        if doresample
            dataPast = dataPastOrig;
            tiFutureDown = 1:length(dataFuture); 
            tiFutureUp = linspace(1, length(dataFuture), samplerat*length(dataFuture));
            dataFuture = interp1(tiFutureDown, dataFuture, tiFutureUp, 'nearest');
            %dataFuture = filter(filtwts,1,dataFuture,filtinit);
            dataFuture = dataFuture';
        end
        % smooth out the kink at past/future interface
        [dataPast, filtinit2] = filter(filtwts,1,dataPast, filtinit);
        dataFuture = filter(filtwts,1,dataFuture,filtinit2);

        % Step 5: Phase and Frequency Estimation
        % Using the Hilbert transform, estimate the current instantaneous
        % phase and frequency, and use this to calculate when the next
        % PhaseOfInterest will likely occur. 
        [t2nextStim,i2nextStim, phEst(tind),frEst(tind)] = blockPDS(...
            dataPast, dataFuture, SamplingFreq, PhaseOfInterest, ...
            filtdelay/SamplingFreq, ... FIR filter imposes this delay (s)
            loco, hico);
        i2nextStim = i2nextStim - filtdelay; % in terms of current time

        % Step 6: Send a stimulus pulse when appropriate 
        i2nextStim_prev = i2nextStim_prev - packetLength; % one packet passed
        if i2nextStim_prev < 0
            toStim(tind+i2nextStim_prev) = true;
        end
        i2nextStim_prev = i2nextStim;

        %{
        % debugging erroneous predictions 
        if tind > filtdelay
            phErrNow = phEst(tind) - phAll(tind-filtdelay);
            phErrNow = abs(phErrNow);
            phErrNow = mod(phErrNow, 2*pi);
            if phErrNow > pi
                phErrNow = 2*pi - phErrNow;
            end
            if norm(dataFuture) > 10*norm(dataPast)
                % 
            end
            if phErrNow > pi/4 % set threshold to begin debugging 
                % do debugging
                if (abs(phEst(tind)-pi/2) < .1) || (abs(phEst(tind)+pi/2) < .1)
                    %keyboard;
                end
            end
        end
        %}
    end
    dur(tind) = toc;
end

% re-align timing; correct for filter delay 
phEst = [phEst(filtdelay:end), nan(1,filtdelay-1)]; 
frEst = [frEst(filtdelay:end), nan(1,filtdelay-1)]; 
W = [W(filtdelay:end,:); nan(filtdelay-1,length(w))];
R = [R(filtdelay:end,:); nan(filtdelay-1,length(r))];
filtdelay = ceil(filtdelay/2);
%toStim = [toStim(filtdelay:end), false(1,filtdelay-1)];
dataOneChannelFilt = [dataOneChannelFilt(filtdelay:end), zeros(1,filtdelay-1)];

%% Part C: Evaluate Real-Time results 
% Compare simulated real-time output with offline-computed ground truth 

% limit signals to central 80% to avoid hilbert edge effects 
t1 = floor(.1*length(t)); t2 = ceil(.9*length(t));
for var = [...
        "t", "dataOneChannel", "dataOneChannelFilt2", "phAll", "frAll", ...
        "dataOneChannelWithArtifact", "isArt", "StimTrainRec", ...
        "dataOneChannelFilt", "phEst", "frEst", "toStim"]
    eval(var+" = "+var+"(t1:t2);");
end
W = W(t1:t2,:); R = R(t1:t2,:);

if showplots

% plot artifact removal 
figure; 
ax(1) = subplot(211); 
plot(t, dataOneChannelWithArtifact, 'k'); 
grid on; hold on; 
plot(t, dataOneChannel, 'b'); 
hold on; plot(t, dataOneChannel.*(isArt), '--r');
if artDur > 0
    title('Artifact Removal'); 
else
    title('Artifact Removal Disabled');
end
ylabel(channelName);
ax(2) = subplot(212); 
plot(t, StimTrainRec); 
ylabel('ainp1');
grid on; linkaxes(ax, 'x'); 

% compare instantaneous phase, frequency: estimated vs actual
%phAll2 = sin(phAll); phEst2 = sin(phEst);
phErr = phEst - phAll; % [-2*pi -> 2*pi];
%phErr = mod(phErr + 2*pi, 2*pi); % [0 -> 2*pi];
%phErr = phErr - 2*pi*(phErr >= pi); % [-pi -> pi]
frErr = frEst - frAll;
figure; sgtitle(channelName);
subplot(2,2,1); plot(phAll, phEst, '.'); 
grid on; title('Phase Accuracy'); 
xlabel('Offline Calc. Phase (rad)'); ylabel('RealTime Pred. Phase (rad)'); 
xlim([-pi pi]); ylim([-pi pi]);
subplot(2,2,3); polarhistogram(phErr); title('Phase Error (RealTime-Offline)');
subplot(2,2,2); plot(frAll, frEst, '.'); 
xlim([loco hico]); ylim([loco hico]);
grid on; title('Frequency Accuracy'); 
xlabel('Offline Calc. Freq. (Hz)'); ylabel('RealTime Pred. Freq. (Hz)'); 
subplot(2,2,4); histogram(frErr); title('Freq. Error (RealTime-Offline)');
ylabel('Count'); xlabel('Freq. Error (Hz)'); grid on;

% show time-series of intended stim 
figure; 
plot(t, dataOneChannelFilt2); hold on; grid on; 
stem(t(toStim), dataOneChannelFilt2(toStim));
title('Stimulus Timing');
legend('Data', 'Stim', 'Location','westoutside')
ylabel(channelName); xlabel('time');

% show true phase of intended stim
figure; sgtitle({channelName, ...
    ['Goal = ',num2str(PhaseOfInterest*180/pi),' degrees']})
subplot(2,2,1); polarhistogram(phAll(toStim), 18); 
title('Actual Phase of Stim'); 
subtitle(['RMSE ',num2str(rms(phAll(toStim)-PhaseOfInterest, 'omitnan'))])
subplot(2,2,2); polarhistogram(phEst(toStim), 18);
title('Online Est. Phase of Stim');
subtitle(['RMSE ',num2str(rms(phEst(toStim)-PhaseOfInterest, 'omitnan'))])
subplot(2,2,3); polarhistogram(phAll(StimTrainRec), 18);
title('Actual Phase of Recorded Stim');
subtitle(['RMSE ',num2str(rms(phAll(StimTrainRec)-PhaseOfInterest, 'omitnan'))])

% show AR mdl changes
W = W(~isnan(phEst),:);
W = [w0; W];
R = R(~isnan(phEst),:);
R = [r0; R];
if stepsize > 0
    figure; sgtitle(channelName);
    subplot(1,2,1); imagesc(W); colorbar; 
    title('Dynamic AR model update');
    xlabel('tap'); ylabel('iteration');
    subplot(1,2,2); imagesc(abs(R)); colorbar;
    title('System Pole Magnitudes');
    xlabel('pole'); ylabel('iteration');
end

end

% process function outputs 
phStim = phAll(toStim);
phAll = phAll(~isnan(phEst)); phEst = phEst(~isnan(phEst));
frAll = frAll(~isnan(frEst)); frEst = frEst(~isnan(frEst));

%% helper(s) 
% gradient descent update of AR model weights 


    % Update based on gradient wrt weights; no guarantee of stability.
    function w = updateWts(w, x, y)
        ypred = w*x;
        E = y-ypred; del = x*E;
        if donorm
            del = del./(x'*x + eps);
        end
        w = w + stepsize*del';
    end


    % Update stably based on gradient wrt roots of c.p. which are
    % forced to have magnitude less than 1.
    function [w, r, rAmpArg, rAng] = updateWtsStbl(w, r, rAmpArg, rAng, x, y)

        % grad err wrt wts
        ypred = w*x;
        E = y-ypred; dEdw = x*E;

        % grad err wrt rts
        N = length(w); K = length(r);
        dwdr = zeros(N,K);
        dwdr(N,:) = 1; 
        for nn = 1:(N-1)
            n = N-nn; s = 1-2*mod(nn,2);
            for k = 1:K
                rr = r([(1:(k-1)),((k+1):end)]);
                dwdr(n,k) = s * SumAllProds(rr, nn);
            end
        end
        dEdr = dEdw' * dwdr;

        % grad err wrt mag, phase
        drdAmp = exp(1i*rAng - rAmpArg)./((1+exp(-rAmpArg)).^2);
        drdAng = 1i*exp(1i*rAng)./(1+exp(-rAmpArg));
        dEdAmp = dEdr .* conj(drdAmp); dEdAng = dEdr .* conj(drdAng);
        dEdAmp = real(dEdAmp); dEdAng = real(dEdAng);

        % iterate grad. desc.
        if donorm
            dEdAmp = -dEdAmp./(x'*x + eps);
            dEdAng = -dEdAng./(x'*x + eps);
        end
        rAmpArg = rAmpArg + stepsize*dEdAmp;
        rAng = rAng + stepsize*dEdAng;
        r = exp(1i*rAng) ./ (1+exp(-rAmpArg));

        % convert rts to wts
        for nn = 0:(N-1)
            n = N-nn; s = 1-2*mod(nn,2);
            w(n) = s * SumAllProds(r, nn+1);
        end
        if rms(imag(w)) > .05*rms(abs(w))
            % this shouldn't happen, because weights should be real
            keyboard
        end
        w = real(w);
    end

%{
    function S = SumAllProds(vals, ord)
        if ord == 1
            S = sum(vals); % base case
        else 
            S = 0;
            for k = 2:length(vals)
                S = S + vals(k-1)*SumAllProds(vals(k:end), ord-1);
            end
        end
    end
%}

function S = SumAllProds(vals, n)
    % Fast computation of sum of all products of combinations of n elements
    % Input:
    %   vals: vector of numeric values
    %   n: number of elements in each combination (positive integer)
    % Output:
    %   S: sum of products of all n-element combinations (no repetitions)

    if n == 0
        S = 1; % Convention: product of 0 elements is 1
        return;
    elseif n > numel(vals)
        S = 0;
        return;
    end

    if n == 1
        S = sum(vals);
        return;
    end

    C = nchoosek(1:numel(vals), n);       % All index combinations
    P = prod(vals(C), 2);                 % Product along each row
    S = sum(P);                           % Sum of all products
end

end
