function [FF, t] = AnimateHand(deviceData)

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

FF = nan(5,5,3,length(t));
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
FF = [WW; PP; FF];

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

end