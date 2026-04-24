function [dtaPred, figout] = AdaptKalmanAuton(noiseRef, dtaTbl, Am, KA, A0, stepsize, Q, N, w0, ...
    nLMS, dLMS, showfit)
% 
% Adaptive LTV sys ID plus artifact removal with Kalman+LMS filters. 
% The least mean squares (LMS) filter adaptively estimates the artifact and
% subtracts it. The LMS output amplitude also provides an estimate of the 
% measurement noise covariance (R) used by the Kalman filter. The Kalman
% filter makes an estimate of the state using matrix A, makes an adaptive
% estimate of the process noise, and smoothly weighs between state estimate
% and measurement-based estimate of signal value. The state matrix itself
% can also be identified adaptively if nonzero/nonempty parameters Am and
% KA are provided to define the adaptive update rule. Current version only
% supports full state availability, i.e. system output = full state (no
% latent states). 
% 
% Inputs: 
%   noiseRef: reference signal for start of artifacts 
%   dtaTbl: input noisy data as timetable 
%   Am: Hurwitz stabilizing matrix for adaptive state matrix ID; 
%       leave empty for no adaptive system ID
%   KA: gain matrix for adaptive state matrix update; 
%       leave empty for no adaptive system ID
%   A0: initial estimate for state matrix 
%   stepsize: learning rate for LMS filter 
%   Q: process noise covariance estimate for Kalman filter
%   N: LMS filter length (number of taps); should be at least as long as
%      expected artifact duration (in samples) 
%   w0: initial LMS filter weights; default = all zero 
%   nLMS: flag to normalize LMS step update; default = true 
%   dLMS: flag to minimize noise derivative instead of noise value; default
%         = false
%   showfit: flag to produce figures and display progress; default = false
% 
% Outputs: 
%   dtaPred: output data with artifact removed 
%   figout: figure handles 
% 

%% arg handling 

% handle incomplete args
if nargin < 12
    showfit = false;
    if nargin < 11
        dLMS = false;
        if nargin < 10
            nLMS = true;
        end
    end
end

%N = size(RR,3);

% handle empty args 
if isempty(Am)
    Am = -eye(width(dtaTbl));
elseif numel(Am) == 1
    Am = Am*eye(width(dtaTbl));
end
if isempty(KA)
    KA = lyap(Am',eye(size(Am)));
elseif numel(KA) == 1
    %KA = lyap(Am',KA*eye(size(Am)));
    KA = KA*eye(size(Am));
end
if isempty(A0)
    A0 = eye(width(dtaTbl)); % discrete
end
if isempty(w0)
    w0 = zeros(N,width(dtaTbl));
end

Th = seconds(dtaTbl.Properties.TimeStep);
if isnan(Th)
    Th = median(seconds(diff(dtaTbl.Time)));
end

%% check for stability 

% check Am Hurwitz 
lambda = eig(Am); 
if sum(lambda >= 0)
    warning('Am should be Hurwitz.')
end
% check KA PDS 
lambda = eig(KA);
if sum(lambda <= 0) | ~isequal(KA,KA')
    warning('KA should be PDS.')
end
% check lyap deriv NDS
dLyap = Am'*KA + KA*Am; 
lambda = eig(-dLyap);
if sum(lambda <= 0) | ~isequal(dLyap,dLyap')
    warning('KA produces unstable system.')
end

%% main run 

% raw vars setup
Y = table2array(dtaTbl)';
Xest = nan(size(Y));
noiseRef = [mode(noiseRef)*ones(1,N-1), noiseRef']; % pad

% AID / Kalman setup
y = Y(:,1); 
AmInv = Am^-1;
Ahat_d = A0;
xest = y; 
Xest(:,1) = xest;
A_d = expm(Am*Th);
Beta_d = Ahat_d - A_d;
M = Th*AmInv*(A_d-eye(size(A_d)));
n = N;
%Q = cov(Y(:,1:N)'); 
P = Q; R = Q; % default/starting assumption

% LMS setup 
w = w0;
Gprev = zeros(1,N);
Eprev = zeros(1,height(Y));
LMSdata = nan(size(dtaTbl));

% track observer noise - NOT USED
ObsNoise = nan(N, width(dtaTbl), 256);
    % dim 1: tap
    % dim 2: data channels 
    % dim 3: # examples to track 
ObsNoiseP = 1; % next index of dim 3

figout = [];
if showfit 
    figout = figure('Units','Normalized','Position',[.1 .1 .8 .8]); 
    subplot(2,2,1); 
    imgL = imagesc(A0); colorbar; 
    title('Disc A LSQ fit');
    cmin = min(A0(:)); cmax = max(A0(:)); cdiff = cmax-cmin; 
    if cdiff == 0
        cdiff = .01;
    end
    cmin = cmin - .05*cdiff; cmax = cmax + .05*cdiff;
    imgL.Parent.CLim = [cmin, cmax];
    subplot(2,2,2); 
    img = imagesc(A0, [cmin,cmax]); colorbar;
    ttxt = title('Disc A Adapt: 0%');
    subplot(2,2,3); 
    wplot = stem(w); grid on; 
    title('LMS online weights'); xlabel('tap'); ylabel('weight');
    subplot(2,2,4); 
    LMSdataplot = plot(LMSdata); grid on;
    title('LMS noise estimate'); xlabel('timepoint'); ylabel('signal');
    drawnow;
    prog_curr = 0; prog_prev = 0; prog_updTick = 1; % percent progress
end

for t = 2:width(Y)

    if noiseRef(t+N-1)
        n = 1;
    else
        n = n+1;
    end
    y = Y(:,t);

    % LMS est observer noise 
    Gidx = noiseRef((1:N)+t-1);
    if sum(abs(Gidx))
        %keyboard
    end
    noiseLMS = Gidx*w;
    EE = y' - noiseLMS;
    if ~dLMS
        dw = Gidx'*EE;
        if nLMS
            dw = dw./(Gidx*Gidx' + eps);
        end
    else
        dG = Gidx-Gprev;
        dw = dG' * (EE-Eprev);
        if nLMS
            dw = dw./(dG*dG' + eps);
        end
    end
    w = w + stepsize*dw;
    Gprev = Gidx; Eprev = EE;
    y = EE'; % subtract LMS estimate of noise 
    ObsNoise(:,:,ObsNoiseP) = w; 
    if ObsNoiseP < size(ObsNoise,3)
        ObsNoiseP = ObsNoiseP + 1;
    else
        ObsNoiseP = 1; % begin overwriting beginning
    end
    LMSdata(t,:) = noiseLMS;

    % adapt update state equation 
    dAhat = -KA*(xest-y)*y';
    delBeta_d = M*dAhat;
    if n > N
        Beta_d = Beta_d + delBeta_d; 
        % update A matrix only when no noise. 
        % CONSIDER: replace this crude rule with a Kalman-like gain based on noise?
    end
    Ahat_d = Beta_d + A_d;

    % observer noise 
    %{
    RR = var(ObsNoise,[],3, 'omitnan'); % variance, tap x channel
    RRnow = Gidx*RR; % variance by channel
    if ~sum(isnan(RRnow(:)))
    R = diag(RRnow);
        % R is observer noise, which depends on the noiseRef signal. 
        % As R vanishes, K becomes I, and xpos becomes y (only observation
        % is used, not state transition estimate). 
    end
    %}
    R = diag(noiseLMS.^2);

    % Kalman predict 
    Ppri = Ahat_d*P*Ahat_d' + Q;
    xpri = Ahat_d*xest;

    % Kalman update     
    K = Ppri * (Ppri + R)^-1;
    xpos = xpri + K*(y-xpri);
    Ppos = (eye(size(K)) - K)*Ppri;

    % Kalman forward project and store 
    P = Ppos;
    xest = xpos;
    Xest(:,t) = xest;

    if showfit
        prog_curr = 100*t/width(Y); 
        if prog_curr - prog_prev > prog_updTick
            prog_prev = prog_curr;
            [cmin,chg_min] = min([cmin, min(Ahat_d(:))]); chg_min = chg_min-1;
            [cmax,chg_max] = max([cmax, max(Ahat_d(:))]); chg_max = chg_max-1;
            img.CData = Ahat_d;
            if chg_max | chg_min
                img.Parent.CLim = [cmin, cmax];
                imgL.Parent.CLim = [cmin, cmax];
            end
            ttxt.String = ['Disc A Estimate: ',num2str(prog_curr,3),'%'];
            for ch = 1:width(dtaTbl)
                wplot(ch).YData = w(:,ch);
                LMSdataplot(ch).YData = LMSdata(:,ch);
            end
            pause(1e-9); drawnow; pause(1e-9); 
        end
    end

end

dtaPred = dtaTbl; 
dtaPred{:,:} = Xest';

end