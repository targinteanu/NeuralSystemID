%% load files 
load('/Users/torenarginteanu/Desktop/Data_PD/PD26N003/baseline_1.mat', 'deviceData1');
NS = openNSx('/Users/torenarginteanu/Desktop/Data_PD/PD26N003/NSdata/ECoG PDS/datafile002.ns2','uV');
NEV = openNEV('/Users/torenarginteanu/Desktop/Data_PD/PD26N003/NSdata/ECoG PDS/datafile002.nev');
NS = ns2timetable(NS); NEV = nev2table(NEV);
t0 = NS.Time(350000); % avoid stim
NS = NS(NS.Time > t0,:); NEV = NEV(NEV.Time > t0,:);
dt = 0.001; % 1/Fs

%% process hand grasping 
[FF, t1, Tthetaphi, Fl] = AnimateHand(deviceData1);
z1 = interp1(t1, deviceData1.f(2).dEnds.z, t1(1):dt:t1(end), "linear","extrap");
z1 = z1.*(z1>0);
Th = Tthetaphi.Variables;
Th = interp1(seconds(Tthetaphi.Time),Th, t1(1):dt:t1(end), "linear","extrap");

%% extract serial data 
NEVsrl = arrayfun(@(lbl) sscanf(lbl, 'SerialDigitalIO: %f'), NEV.EventLabels);
tNEV = NEV.Time; 
tNEV = seconds(tNEV - tNEV(1));
tNEV = tNEV(NEVsrl < 253); NEVsrl = NEVsrl(NEVsrl < 253); % ignore stop/go prompts
zNEV = interp1(tNEV,NEVsrl, tNEV(1):dt:tNEV(end), "linear", "extrap");
tNEV = tNEV(1):dt:tNEV(end);

%% align grasp/Ephys 

[R,L] = xcorr(z1, zNEV);
[r,l] = max(R); l = L(l);
if l>0
    z1 = z1((l+1):end);
    Th = Th((l+1):end,:);
    t = tNEV;
elseif l<0
    zNEV = zNEV((1-l):end);
    t = tNEV((1-l):end);
end

dur = min(height(Th), length(zNEV));
Th = Th(1:dur, :);
t = t(1:dur);

X = NS.Variables;
X = interp1(seconds(NS.Time - t0), X, t, "linear","extrap");
X = X(:,1:63);

%% select relevant channels 
% should be based on anatomy instead of signal! 
B = X\Th; b = rms(B,2);
bthresh = quantile(b, .75);
X = X(:, b>bthresh);

%% inspect features 

vel = movmean(abs(diff(rms(Th,2))), ceil(1/dt));
[~,movthresh] = midcross(vel); movon = vel > movthresh;

SS = cell(4,width(X));
for ch = 1:width(X)
    [S,Sf,St] = spectrogram(X(:,ch), ceil(1/dt), ceil(1/dt)-1, [], 1/dt);
    S = abs(S);
    movon_ch = movon((St(1)/dt):(St(end)/dt));
    Son = S(:,movon_ch); Soff = S(:,~movon_ch);
    SS{1,ch} = mean(Son, 2); SS{2,ch} = std(Son, [],2).^2;
    SS{3,ch} = mean(Soff,2); SS{4,ch} = std(Soff,[],2).^2;
end

Son_avg  = mean(cell2mat(SS(1,:)), 2);
Son_var  = mean(cell2mat(SS(2,:)), 2);
Soff_avg = mean(cell2mat(SS(3,:)), 2);
Soff_var = mean(cell2mat(SS(4,:)), 2);

figure; 
errorbar(Sf, Son_avg,  -sqrt(Son_var),  sqrt(Son_var));
hold on; grid on; 
errorbar(Sf, Soff_avg, -sqrt(Soff_var), sqrt(Soff_var));
set(gca(), 'YScale', 'log');
ylabel('freq (hz)'); xlabel('power'); title('power spectrum');
legend('Movement On', 'Movement Off')

%% process features 
BPF1 = fir1(1023,[2,12]*2*dt);
BPF2 = fir1(1023,[15,25]*2*dt);
X = [filtfilt(BPF1,1,X), filtfilt(BPF2,1,X)];

%% decoding (standard)
trn = ceil(height(X)/2);

% pick a few fingers to decode 
fnames = Tthetaphi.Properties.VariableNames;
Y = Th(:, ...
    ~contains(fnames, 'j4') & ...
    (contains(fnames, 'd2') | contains(fnames, 'd3')));

% compute output sigmoid nonlinearity 
Y = tan(Y-median(Y));

% complete decoding 
X_aug = [X, ones(height(X),1)];
X1 = X_aug(1:trn,:); X2 = X_aug((trn+1):end,:);
Y1 = Y(1:trn,:); Y2 = Y((trn+1):end,:);
b = X1\Y1;
Y1e = X1*b; Y2e = X2*b;
figure; plot(Y(:,2)); hold on; grid on; plot([Y1e(:,2);Y2e(:,2)]); plot(Y1e(:,2));

disp('Train Acc:')
disp(diag(corr(Y1, Y1e)))
disp('Testing Acc:')
disp(diag(corr(Y2, Y2e)))

%% ID state space system 

Xa = X(1:(end-1),:); Xb = X(2:end,:);
Ya = Y(1:(end-1),:); Yb = Y(2:end,:);
X1a = Xa(1:trn,:); X2a = Xa((trn+1):end,:);
X1b = Xb(1:trn,:); X2b = Xb((trn+1):end,:);
Y1a = Ya(1:trn,:); Y2a = Ya((trn+1):end,:);
Y1b = Yb(1:trn,:); Y2b = Yb((trn+1):end,:);
X1a_aug = [X1a,Y1a]; X1b_aug = [X1b,Y1b]; 
X2a_aug = [X2a,Y2a]; X2b_aug = [X2b,Y2b]; 

AT = X1a_aug\X1b_aug;
C = eye(size(AT)); C = C(1:width(X),:);
sysd = ss(AT',[],C,0,dt);
sysc = d2c(sysd);
O = obsv(sysc);
disp("Full system observability is rank "+rank(O)+" out of "+width(O));

%% build LQE 

% reduced order
Auu = sysc.A((width(X)+1):end, (width(X)+1):end);
Aau = sysc.A(1:(width(X)), (width(X)+1):end);
Aua = sysc.A((width(X)+1):end, 1:width(X));
Aaa = sysc.A(1:(width(X)), 1:width(X));
O = obsv(Auu, Aau);
disp("Reduced order observability is rank "+rank(O)+" out of "+width(O));
R = eye(size(Aau,1)); % penalty term
L = lqe(Auu,zeros(height(Auu),0),Aau,[],R);

% dxu/dt ~ (Auu-LAau)*xu + (Aua-LAaa)*xa + L*dxa/dt
% w := xu - L*xa
% dw/dt ~ (Auu-LAau)*w + ((Aua-LAaa)+(Auu-L*Aau)*L)*xa
obsc = ss(...
    Auu - L*Aau, ...
    Aua - L*Aaa + (Auu-L*Aau)*L, ...
    eye(height(Auu)), 0);
obsd = c2d(obsc, dt);

% decoding (LQE-based) 
W1b = zeros(width(Y),height(X1b)); W1b(:,1) = Y1b(1,:)'   - L*X1b(1,:)';
W2b = zeros(width(Y),height(X2b)); W2b(:,1) = Y1b(end,:)' - L*X1b(end,:)';
for ti = 2:width(W1b)
    W1b(:,ti) = obsd.A*W1b(:,ti-1) + obsd.B*X1b(ti-1,:)';
end
for ti = 2:width(W2b)
    W2b(:,ti) = obsd.A*W2b(:,ti-1) + obsd.B*X2b(ti-1,:)';
end
Y1be = W1b + L*X1b'; Y1be = Y1be';
Y2be = W2b + L*X2b'; Y2be = Y2be';

figure; plot(Yb(:,2)); hold on; grid on; 
plot([Y1be(:,2);Y2be(:,2)]); plot(Y1be(:,2));

disp('Train Acc:')
disp(diag(corr(Y1b, Y1be)))
disp('Testing Acc:')
disp(diag(corr(Y2b, Y2be)))