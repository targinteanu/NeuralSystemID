function [FF, t] = AnimateHand(deviceData)
%% translate struct to matrix of raw finger pos

t = deviceData.time; dt = [0,diff(t)];

W = [deviceData.wristPosition.x; 
     deviceData.wristPosition.y;
     deviceData.wristPosition.z];
P = [deviceData.palmPosition.x; 
     deviceData.palmPosition.y;
     deviceData.palmPosition.z];

F = deviceData.f;
F = struct2cell(F);
F = squeeze(F);
F = F([6,2,3,4,5],:);

FF = nan(5,5,3,length(t)); % joint, finger, XYZ, t
for f = 1:5
    for j = 1:5
        FF(j,f,1,:) = F{j,f}.x;
        FF(j,f,2,:) = F{j,f}.y;
        FF(j,f,3,:) = F{j,f}.z;
    end
end

WW = nan([1,1,size(W)]); WW(:) = W(:);
WW = repmat(WW, [1,5,1,1]);
PP = nan([1,1,size(P)]); PP(:) = P(:);
PP = repmat(PP, [1,5,1,1]);
%FF = [WW; PP; FF];
FF = [WW; FF] - PP; % place palm at origin

%% outlier detection 
% throw out data where hand jumps impossibly 

%% animate 

myfig = figure; 

for f = 1:5
    Fp(f) = plot3(FF(:,f,1,1), FF(:,f,2,1), FF(:,f,3,1), '-o');
    hold on; grid on;
end
x = FF(:,:,1,:); x = x(:);
y = FF(:,:,2,:); y = y(:);
z = FF(:,:,3,:); z = z(:);
xlim([min(x), max(x)]);
ylim([min(y), max(y)]);
zlim([min(z), max(z)]);
xlabel('x'); ylabel('y'); zlabel('z');
mytitle = title(['Time: 0 of ',num2str(t(end)),' seconds.']);

dtskip = 0;
for ti = 2:size(FF,4)
    dti = dt(ti) - dtskip; 
    dti = min(1,dti); dti = max(0,dti);
    pause(dti);
    startTic = tic;
    if isvalid(myfig)
        mytitle.String = ['Time: ', num2str(t(ti)), ' of ', num2str(t(end)), ' seconds.'];
        for f = 1:5
            Fp(f).XData = FF(:,f,1,ti);
            Fp(f).YData = FF(:,f,2,ti);
            Fp(f).ZData = FF(:,f,3,ti);
        end
        drawnow; 
    else
        warning('Stopping animation due to closed figure.')
        break
    end
    dtskip = toc(startTic);
end

%% translate raw finger pos to useful measurement 

Fd = diff(FF, 1); % joint displacement 
Fl = squeeze(sqrt(sum(Fd.^2,3))); % joint length
Fl = median(Fl,3);

% digital flexion/extension; thenar abduction/adduction
Ftheta = nan(size(Fd,1)-1, size(Fd,2), size(Fd,3));
for j = 2:size(Fd,1)
    for f = 1:size(Fd,2)
        for ti = 1:size(Fd,4)
            dp = squeeze(Fd(j,f,:,ti))' * squeeze(Fd(j-1,f,:,ti));
            costheta = dp / (Fl(j,f) * Fl(j-1,f));
            costheta = max(-1, costheta); costheta = min(1, costheta);
            Ftheta(j-1,f,ti) = acos(costheta); 
        end
    end
end
Ftheta(2,1,:) = nan; % thumb lacks final joint 
Ftheta(1,:,:) = nan; % no real movement m wrt wrist 
Ftheta(end,2:end,:) = nan; % final joint not independent

% wrist flexion/extension, pron/supination, deviation
Pd = [Fd(2,3,:,:); FF(3,2,:,:)-FF(3,3,:,:)]; % middle mcp and index-middle knuckle disp
Pd = squeeze(Pd);
Ptheta = nan(3,size(Pd,3));
for ti = 1:size(Pd,3)
    xp = cross(squeeze(Pd(1,:,ti)), squeeze(Pd(2,:,ti))); % palm normal
    [Ptheta(1,ti), Ptheta(2,ti)] = cart2sph(xp(1),xp(2),xp(3));
    Ptheta(3,ti) = cart2sph(Pd(1,1,ti),Pd(1,2,ti),Pd(1,3,ti));
end

FD = FF(end,:,:,:) - FF(1,:,:,:); % finger displacement 
FD = squeeze(FD);
FL = squeeze(sqrt(sum(FD.^2,3))); % finger length
FL = median(FL,2);
FD3 = squeeze(FD(3,:,:)); FD = FD([1,2,4,5],:,:); 
FL3 = FL(3); FL = FL([1,2,4,5]);

% digital absuction/adduction; thenar flexion/extension
Fphi = nan(size(FD,1),size(FD,3));
for f = 1:size(FD,1)
    for ti = 1:size(FD,3)
        dp = squeeze(FD(f,:,ti)) * FD3(:,ti);
        costheta = dp / (FL(f) * FL3); % calculate cosine of angle
        costheta = max(-1, min(1, costheta)); % clamp value for acos
        Fphi(f, ti) = acos(costheta); 
    end
end

%% plot 
Fthetaphi = [reshape(Ftheta, [size(Ftheta,1)*size(Ftheta,2),size(Ftheta,3)]); Fphi; Ptheta];
Fthetaphi = Fthetaphi';
Fthetaphi = Fthetaphi(:, any(~isnan(Fthetaphi))); 
Fthetaphi = Fthetaphi(dt>0,:); tthetaphi = t(dt>0);
Ts = min(dt(dt>0)); % resample time 
treg = tthetaphi(1):Ts:tthetaphi(end);
Freg = interp1(tthetaphi,Fthetaphi,treg);
figure; plot(treg, Freg); grid on; xlabel('time (s)'); ylabel('angle (rad)');
figure; pwelch(Freg,[],[],[],1/Ts,'power'); 

end