function dtaPred = AdaptKalmanAuton(noiseRef, dtaTbl, Am, KA, A0, stepsize, Q, N, w0, ...
    nLMS, dLMS, showfit)
% Adaptive LTV sys ID plus artifact removal with Kalman filter 
% TO DO: replace RR with N and have QQ calculated using the output of a LMS
% filter 

% is Am actually being used??

%% arg handling 

% handle incomplete args 
if nargin < 11
    dLMS = false;
    if nargin < 10
        nLMS = true;
        if nargin < 12
            showfit = false;
        end
    end
end

%N = size(RR,3);

% handle empty args 
if isempty(Am)
    Am = -eye(width(dtaTbl));
end
if isempty(KA)
    KA = lyap(Am',eye(size(Am)));
end
if isempty(A0)
    A0 = eye(width(dtaTbl));
end
if isempty(w0)
    w0 = zeros(N,width(dtaTbl));
end

Th = seconds(dtaTbl.Properties.TimeStep);
if isnan(Th)
    Th = mean(seconds(diff(dtaTbl.Time)));
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

% Kalman setup
y = Y(:,1); 
A = A0;
xest = y; 
Xest(:,1) = xest;
[Ad,Bd0] = c2d(Am,eye(size(Am)),Th);
Ac0 = d2c(A,zeros(size(A,1),1),Th);
Ac = Ac0;
n = N;
%Q = cov(Y(:,1:N)'); 
P = Q; R = Q; % default/starting assumption

% LMS setup 
w = w0;
Gprev = zeros(1,N);
Eprev = zeros(1,height(Y));
LMSdata = nan(size(dtaTbl));

% track observer noise 
ObsNoise = nan(N, width(dtaTbl), 256);
    % dim 1: tap
    % dim 2: data channels 
    % dim 3: # examples to track 
ObsNoiseP = 1; % next index of dim 3

if showfit 
    figure('Units','Normalized','Position',[.1 .1 .8 .8]); 
    subplot(2,2,1); 
    imgL = imagesc(Ac); colorbar; 
    title('Cont A LSQ fit');
    cmin = min(Ac(:)); cmax = max(Ac(:)); cdiff = cmax-cmin; 
    if cdiff == 0
        cdiff = .01;
    end
    cmin = cmin - .05*cdiff; cmax = cmax + .05*cdiff;
    imgL.Parent.CLim = [cmin, cmax];
    %Ac = zeros(size(Ac)); % reset 
    subplot(2,2,2); 
    img = imagesc(Ac, [cmin,cmax]); colorbar;
    ttxt = title('Cont A Adapt: 0%');
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
    dA = -KA*(xest-y)*y';
    if n > N
        Ac = Ac + dA*Th; 
        % update A matrix only when no noise. 
        % CONSIDER: replace this crude rule with a Kalman-like gain based on noise?
    end
    %[Ad,Bd] = c2d(Am,Ac-Am,Th);
    Bd = Bd0*(Ac-Am);
    Ad2 = expm(Ac*Th);

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
    %Ppri = Ad*P*Ad' + Bd*(y*y')*Bd' + Q;
    %xpri = Ad*xest + Bd*y;
    Ppri = Ad2*P*Ad2' + Q;
    xpri = Ad2*xest;

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
            [cmin,chg_min] = min([cmin, min(Ac(:))]); chg_min = chg_min-1;
            [cmax,chg_max] = max([cmax, max(Ac(:))]); chg_max = chg_max-1;
            img.CData = Ac;
            if chg_max | chg_min
                img.Parent.CLim = [cmin, cmax];
                imgL.Parent.CLim = [cmin, cmax];
            end
            ttxt.String = ['Cont A Estimate: ',num2str(prog_curr,3),'%'];
            for ch = 1:width(dtaTbl)
                wplot(ch).YData = w(:,ch);
                LMSdataplot(ch).YData = LMSdata(:,ch);
            end
            drawnow;
        end
    end

end

dtaPred = dtaTbl; 
dtaPred{:,:} = Xest';

end