function dtaPred = AdaptKalmanAuton(noiseRef, dtaTbl, Am, KA, A0, mu, QQ, w0, showfit)
% Adaptive LTV sys ID plus artifact removal with Kalman filter 
% TO DO: replace QQ with N and have QQ calculated using the output of a LMS
% filter 

%% arg handling 

% handle incomplete args 
if nargin < 9
    showfit = false;
end

N = size(QQ,3);

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

Y = table2array(dtaTbl)';
Xest = nan(size(Y));

y = Y(:,1); 
A = A0;
xest = y; 
Xest(:,1) = xest;
[Ad,Bd0] = c2d(Am,eye(size(Am)),Th);
Ac0 = d2c(A,zeros(size(A,1),1),Th);
Ac = Ac0;
n = N;
Q = cov(Y(:,1:N)'); P = Q;

if showfit 
    figure('Units','Normalized','Position',[.1 .3 .8 .4]); 
    subplot(121); imgL = imagesc(Ac); colorbar; 
    title('Cont A LSQ fit');
    cmin = min(Ac(:)); cmax = max(Ac(:)); cdiff = cmax-cmin; 
    if cdiff == 0
        cdiff = .01;
    end
    cmin = cmin - .05*cdiff; cmax = cmax + .05*cdiff;
    imgL.Parent.CLim = [cmin, cmax];
    Ac = zeros(size(Ac)); % reset 
    subplot(122); img = imagesc(Ac, [cmin,cmax]); colorbar;
    ttxt = title('Cont A Adapt: 0%');
    pause(1e-9);
    prog_curr = 0; prog_prev = 0; prog_updTick = 1; % percent progress
end

for t = 2:width(Y)

    if noiseRef(t)
        n = 1;
    else
        n = n+1;
    end

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

    % predict 
    %Ppri = Ad*P*Ad' + Bd*(y*y')*Bd' + Q;
    %xpri = Ad*xest + Bd*y;
    Ppri = Ad2*P*Ad2' + Q;
    xpri = Ad2*xest;

    % update 
    R = QQ(:,:,min(n,N));
    K = Ppri + (Ppri + R)^-1;
    y = Y(:,t);
    xpos = xpri + K*(y-xpri);
    Ppos = (eye(size(K)) - K)*Ppri;

    % forward project and store 
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
            pause(1e-9);
        end
    end

end

dtaPred = dtaTbl; 
dtaPred{:,:} = Xest';

end