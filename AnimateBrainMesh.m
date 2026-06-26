ax = gca();
View0 = ax.View;

%%
ax.View = View0;
dur = 10; % seconds 
playbackspeed = 1; % relative to real time
FR = 60; % frame rate (FPS)
t = 0:(1/FR):dur; % seconds
azpath = View0(1) + 60*sin(2*pi*t/dur);
elpath = View0(2) + 30*sin(4*pi*t/dur);
pausetime = playbackspeed / FR;

VW = VideoWriter('AnimatedBrainMesh', 'MPEG-4');
VW.Quality = 95;
VW.FrameRate = FR;
open(VW);

% setup for cropping
Vframe = getframe(ax.Parent);
VframeSize = size(Vframe.cdata);
pixH = (floor(VframeSize(2)/2):VframeSize(2));
pixV = 1:VframeSize(1);

for ti = 1:length(t)
    ticStart = tic;
    ax.View = [azpath(ti), elpath(ti)];
    Vframe = getframe(ax.Parent);
    Vframe.cdata = Vframe.cdata(pixV,pixH,:);
    writeVideo(VW,Vframe);
    ticDur = toc(ticStart);
    if ticDur < pausetime
        pause(pausetime-ticDur);
    end
end

close(VW);