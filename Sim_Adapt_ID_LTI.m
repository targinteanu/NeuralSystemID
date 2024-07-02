Asign = [-1, +1, +1,  0,  0,  0,  0;
          0, -1, +1,  0,  0,  0,  0;
          0,  0, -1, -1,  0, -1,  0;
          0,  0,  0, -1, -1, -1,  0;
          0,  0,  0,  0, -1, +1,  0;
          0,  0,  0,  0,  0, -1, -1;
         +1,  0,  0,  0,  0,  0, -1];

Bsign = [+1, 0, 0, 0,  0,  0, 0;
          0, 0, 0, 0, -1,  0, 0;
          0, 0, 0, 0,  0, -1, 0]';

C = [1,0,0,0,0,0,0;
     0,0,0,0,1,0,0];

%% setup true system with random coeffs

maxcoeff = 1; 

A = Asign; B = Bsign;

for r = 1:size(A,1)
    for c = 1:size(A,2)
        A(r,c) = maxcoeff*rand*A(r,c);
    end
end
for r = 1:size(B,1)
    for c = 1:size(B,2)
        B(r,c) = maxcoeff*rand*B(r,c);
    end
end

observability_check = rank(obsv(A,C))
eig(A)

%% simulate dynamics 

t = 0:.001:30; % s 
u = rand([length(t),size(B,2)]) > .99; % random impulse train 
u = double(u); 

sys = ss(A,B,C,zeros(size(C,1),size(B,2)));
[y,t,x] = lsim(sys,u,t);

% init values 
Ai = Asign; Bi = Bsign;
xObsi = zeros(size(x,2),1);
L = lqe(Ai,zeros(size(Ai,1),1),C,0,eye(size(C,1)));
xObs = zeros(size(xObsi,1),length(t)); xObs(:,1) = xObsi; 
Aobs = zeros([size(Ai),length(t)]); Aobs(:,:,1) = Ai;
Bobs = zeros([size(Bi),length(t)]); Bobs(:,:,1) = Bi;
LyaV = zeros(size(t));

for ti = 2:length(t)
    dt = t(ti) - t(ti-1);
    yi = y(ti,:)'; ui = u(ti,:)';
    yObsi = C*xObsi; dely = yi - yObsi;
    delx = x(ti,:)' - xObsi;
    delA = Ai - A; delB = Bi - B;

    dxObs = (Ai-L*C)*xObsi + Bi*ui + L*yi;
    dA = C' * dely * xObsi'; 
    dB = C' * dely * ui';

    Ai = Ai + dA*dt; Bi = Bi + dB*dt; 
    xObsi = xObsi + dxObs*dt; 
    xObs(:,ti) = xObsi; 
    Aobs(:,:,ti) = Ai; Bobs(:,:,ti) = Bi;
    LyaV(ti) = .5*( dely'*dely + trace(delA*delA') + trace(delB*delB') );
end
yObs = C*xObs;
LyaV = LyaV(2:end);

%% display results 
figure; 
colr = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F", "g", "c"];

subplot(3,1,1); 
for v = 1:size(x,2)
    plot(t,x(:,v),   '-', 'Color',colr(v),'LineWidth',1);
    hold on; 
    plot(t,xObs(v,:),':', 'Color',colr(v),'LineWidth',1.25);
end
vy0 = v;
for v = 1:size(y,2)
    plot(t,y(:,v),   '-', 'Color',colr(v+vy0),'LineWidth',2);
    hold on; 
    plot(t,yObs(v,:),':', 'Color',colr(v+vy0),'LineWidth',2.25);
end
grid on; xlabel('t (s)'); 

subplot(3,1,2);
Aerr = A-Aobs; 
Aerr = Aerr.^2; Aerr = sum(sum(Aerr)); Aerr = squeeze(Aerr);
Berr = B-Bobs; 
Berr = Berr.^2; Berr = sum(sum(Berr)); Berr = squeeze(Berr);
plot(t,Aerr); hold on; plot(t,Berr); 
legend('A','B'); ylabel('error'); xlabel('t (s)'); 
grid on;

subplot(3,1,3); 
plot(t(2:end),LyaV); grid on; hold on;
plot(t(3:end), diff(LyaV)./diff(t(2:end)));
legend('V', 'dV/dt');
ylabel('Lyap'); xlabel('t (s)');