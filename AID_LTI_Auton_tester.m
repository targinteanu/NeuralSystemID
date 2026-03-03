%% randomly generate a nearly-unstable system 

% transfer function 
numpoles = 8; 
mysysc = tf(1);
for p = 1:numpoles
    freq = 2*pi*50*rand; 
    damp = rand-0.9; 
    mysysc = mysysc * tf(freq, [1, -2*damp, damp^2 + freq^2]);
end
figure; impulse(mysysc); 

% state space 
Ts = .001; % seconds
mySSc = ss(mysysc);
mySSd = c2d(mySSc, Ts);
Ac = mySSc.A; Ad = mySSd.A;

%% simulate autonomous system 
t = 0:Ts:10; % time vector for simulation
n = size(Ac, 1); % state dimension 
x0 = rand(n, 1); % initial state
X = zeros(size(x0,1), length(t)); X(:,1) = x0; x = x0;
for ti = 2:length(t)
    x = Ad*x; % state update 
    X(:,ti) = x + 0.1*randn(size(x)); % measurement noise
end
plot(t, X); % plot the output response
xlabel('Time (s)');
ylabel('States');
title('Response of the Autonomous System');

%% ID discrete 

% estimate overall Ad using Lsq 
X2T = X(:,2:end)'; X1T = X(:,1:(end-1))';
% Ad * X1 = X2 ; X1T * AdT = X2T
Ad_lsq = (X1T \ X2T)'; % Estimate Ad using least squares

% estimate windowed Ad using Lsq
winlen = round(0.5 / Ts); % samples 
numwin = floor(length(t)/winlen);
Ad_win = zeros([size(Ad_lsq), numwin]);
for w = 1:numwin
    startIdx = (w-1)*winlen + 1;
    endIdx = min(w*winlen, length(t)-1);
    Ad_win(:,:, w) = (X1T(startIdx:endIdx, :) \ X2T(startIdx:endIdx, :))';
end

% AID 
    % define params 
    KA = (1e-4)*eye(n);
    Am = (-1e3)*eye(n);
    Q = -eye(n); 
    P = lyap(Am', Q);
Xtbl = array2timetable(X', "RowTimes",seconds(t));
[~,Xtbl_AID,~,AIDeval,~,Ad_AID] = AID_LTI_auton([],Xtbl,Am,KA,[],true);

% compare A matrices
figure('Units','normalized', 'Position',[.05,.05,.9,.9]);
sgtitle('Discrete A');
subplot(2,3,1); imagesc(Ad); colorbar; title('Actual');
subplot(2,3,4); imagesc(Ad_lsq); colorbar; title('Least Sq. Est.');
subplot(2,3,2); imagesc(mean(Ad_win,3)); colorbar; title('Windowed L.Sq. Est.');
subplot(2,3,5); imagesc(std(Ad_win,[],3)); colorbar; title('Windowed S.D.');
subplot(2,3,3); imagesc(Ad_AID(:,:,end)); colorbar; title('Adaptive Final');
subplot(2,3,6); imagesc(std(Ad_AID,[],3)); colorbar; title('Adaptive S.D.');

%% AID continuous 

y0 = zeros(n+n^2,1);
[ODEt, ODEy] = ode45(@contODEdyn,t,y0);

% unstack 
X_ode = ODEy(1:n,:); bstack_ode = ODEy((n+1):end,:);
Ac_ode = zeros(size(Ac),length(ODEt));
Ac_ode(:) = bstack_ode(:); % "beta"
Ac_ode = Ac_ode + Am; % "Ahat"

% compare A matrices 
figure; sgtitle('Continuous A');
subplot(1,3,1); imagesc(Ac); colorbar; title('Actual');
subplot(1,3,2); imagesc(Ac_ode(:,:,end)); colorbar; title('Final Est');
subplot(1,3,3); imagesc(std(Ac_ode,[],3)); colorbar; title('Adaptive S.D.');


function dy_dt = contODEdyn(ti, yi)

% find true values of x at time ti
[~,i] = min(abs(t-ti));
xtrue = X(:,i);

% unstack y -> xhat, beta
xhat = yi(1:n); bstack = yi((n+1):end);
beta = zeros(n); beta(:) = bstack;

% continuous time dynamics 
dxhat_dt = Am*xhat + beta*xtrue;

% param update
db_dt = -KA*P*(xhat-xtrue)*xtrue';

% restack dx, db -> y
dy_dt = [dxhat_dt; db_dt(:)];

end