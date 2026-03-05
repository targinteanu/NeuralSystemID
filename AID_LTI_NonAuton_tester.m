%% randomly generate a nearly-unstable system 

% transfer function 
numpoles = 8; 
mysysc = tf(1);
for p = 1:numpoles
    freq = 2*pi*50*rand; 
    damp = rand-0.99; 
    mysysc = mysysc * tf(freq, [1, -2*damp, damp^2 + freq^2]);
end
figure; impulse(mysysc); 

% state space 
Ts = .001; % seconds
mySSc = ss(mysysc);
Ac = mySSc.A; 
Bc = (2*rand(size(Ac,1),1)-1) * max(abs(Ac(:)));
[Ad,Bd] = c2d(Ac,Bc, Ts);

% randomize an input signal 
t = 0:Ts:10; % time vector for simulation
u = double(rand(size(t)) > .7); % 1D random impulse train

%% simulate autonomous system 
n = size(Ac, 1); % state dimension 
x0 = rand(n, 1); % initial state
X = zeros(size(x0,1), length(t)); X(:,1) = x0; x = x0;
for ti = 2:length(t)
    x = Ad*x + Bd*u(ti-1); % state update 
    X(:,ti) = x + 0.1*randn(size(x)); % measurement noise
end
figure; clear ax;
ax(1) = subplot(2,1,1);
plot(t, X); grid on; % plot the output response
xlabel('Time (s)');
ylabel('States');
title('Response of the Autonomous System');
ax(2) = subplot(2,1,2); 
plot(t, u); grid on; ylabel('Input'); xlabel('Time (s)');
linkaxes(ax, 'x');

% choose channel(s) to display 
Xdiff = X - mean(X); Xdiff = rms(Xdiff, 2); % diff from mean
[~,chandispind] = min(Xdiff); % most average channel

%% ID discrete 

% estimate overall Ad using Lsq 
X2T = X(:,2:end)'; X1uT = [X(:,1:(end-1)); u(1:(end-1))]';
% Ad * X1 = X2 ; X1T * AdT = X2T
ABd_lsq = (X1uT \ X2T)'; % Estimate Ad using least squares
Ad_lsq = ABd_lsq(:,1:(end-1)); Bd_lsq = ABd_lsq(:,end);

% estimate windowed Ad using Lsq
winlen = round(0.5 / Ts); % samples 
numwin = floor(length(t)/winlen);
Ad_win = zeros([size(Ad_lsq), numwin]); Bd_win = zeros([size(Bd_lsq), numwin]);
for w = 1:numwin
    startIdx = (w-1)*winlen + 1;
    endIdx = min(w*winlen, length(t)-1);
    ABd_win = (X1uT(startIdx:endIdx, :) \ X2T(startIdx:endIdx, :))';
    Ad_win(:,:,w) = ABd_win(:,1:(end-1)); 
    Bd_win(:,:,w) = ABd_win(:,end);
end

% AID 
    % define params 
    KA = (1e-8)*eye(n);
    KB = (1e-4)*eye(n);
    Am = (-1e2)*eye(n);
    %Q = (1e-3)*eye(n); 
    Q = cov(X');
    P = lyap(Am', Q);
Xtbl = array2timetable(X', "RowTimes",seconds(t));
chandispname = Xtbl.Properties.VariableNames{chandispind};
%[~,Xtbl_AID,~,AIDeval,~,Ad_AID] = AID_LTI_auton([],Xtbl,Am,KA*P,[],...
%    chandispname);

% compare A matrices
figure('Units','normalized', 'Position',[.05,.05,.9,.9]);
sgtitle('Discrete A | B');
subplot(2,3,1); imagesc([Ad,Bd]); colorbar; title('Actual');
subplot(2,3,4); imagesc([Ad_lsq,Bd_lsq]); colorbar; title('Least Sq. Est.');
subplot(2,3,2); imagesc(mean([Ad_win,Bd_win],3)); colorbar; title('Windowed L.Sq. Est.');
subplot(2,3,5); imagesc(std([Ad_win,Bd_win],[],3)); colorbar; title('Windowed S.D.');
%subplot(2,3,3); imagesc(Ad_AID(:,:,end)); colorbar; title('Adaptive Final');
%subplot(2,3,6); imagesc(std(Ad_AID,[],3)); colorbar; title('Adaptive S.D.');

%% AID continuous 

y0 = zeros(n+n^2+n,1);
[ODEt, ODEy] = ode45(@(ti,yi) contODEdyn(ti,yi,t,X,u,Am,KA,KB,P), t, y0);
ODEt = ODEt'; ODEy = ODEy';

% unstack 
X_ode = ODEy(1:n,:); bstack_ode = ODEy((n+1):end,:);
Beta_ode = zeros([size([Ac,Bc]),length(ODEt)]);
Beta_ode(:) = bstack_ode(:);
Ac_ode = Beta_ode(:,1:n,:) + Am; % "Ahat"
Bc_ode = Beta_ode(:,(n+1):end,:); % "Bhat"

% compare A matrices 
figure; sgtitle('Continuous A | B');
subplot(1,3,1); imagesc([Ac,Bc]); colorbar; title('Actual');
subplot(1,3,2); imagesc([Ac_ode(:,:,end),Bc_ode(:,:,end)]); colorbar; title('Final Est');
subplot(1,3,3); imagesc(std([Ac_ode,Bc_ode],[],3)); colorbar; title('Adaptive S.D.');

% compare X values 
figure; clear ax;
ax(1) = subplot(2,1,1); 
plot(t, X(chandispind,:)); hold on; grid on; plot(ODEt, X_ode(chandispind,:));
legend('Actual', 'Estimate'); ylabel(chandispname); xlabel('Time (sec)');
x_ode = interp1(ODEt, X_ode(chandispind,:), t, 'linear', 'extrap');
ax(2) = subplot(2,1,2); semilogy(t, (X(chandispind,:)-x_ode).^2); grid on; 
title('Squared Error');
linkaxes(ax, 'x');


function dy_dt = contODEdyn(ti, yi, t, X, u, Am, KA, KB, P)

% find true values of x at time ti
[~,i] = min(abs(t-ti));
xtrue = X(:,i);
n = size(X,1);

% unstack y -> xhat, beta
xhat = yi(1:n); bstack = yi((n+1):end);
beta = zeros(n,n+1); beta(:) = bstack;

% continuous time dynamics 
dxhat_dt = Am*xhat + beta*[xtrue; u(i)];

% param update
db_dt = [-KA*P*(xhat-xtrue)*(xtrue'), -KB*P*(xhat-xtrue)*(u(i))'];

% restack dx, db -> y
dy_dt = [dxhat_dt; db_dt(:)];

end