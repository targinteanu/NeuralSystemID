%% define files of interest
%{
fp = '/Users/torenarginteanu/Desktop/Data_PD/PD25N008/Motor Test';
graspbl = 'Grasping1_125140.mat';
graspstim = 'Grasping1_125033.mat';
%}
%%{
fp = '/Users/torenarginteanu/Desktop/Data_PD/PD26N002/Grasping';
graspbl = 'Grasping2_131426.mat';
graspstim = 'Grasping2_131830.mat';
%}

Fs = 1000; % sample rate of interest

%% 
load(fullfile(fp,graspbl));
graspblData = deviceData2; clear deviceData2
load(fullfile(fp, graspstim));
graspstimData = deviceData2; clear deviceData2

%%
[~,~,T1] = AnimateHand(graspblData);
[~,~,T2] = AnimateHand(graspstimData);

%%
filttrem = buildFIRBPF(Fs,3,7); % parkinsonian tremor 
filtgrasp = buildFIRBPF(Fs,1,3);

selFlEx = contains(T1.Properties.VariableNames, 'j2FlEx'); 

T1r = retime(T1, 'regular', 'linear', 'SampleRate',Fs);
T2r = retime(T2, 'regular', 'linear', 'SampleRate',Fs);
T1i = zeros(size(T1.Time)); T2i = zeros(size(T2.Time));
for i = 1:height(T1)
    [~,T1i(i)] = min(abs(seconds(T1.Time(i)-T1r.Time)));
end
for i = 1:height(T2)
    [~,T2i(i)] = min(abs(seconds(T2.Time(i)-T2r.Time)));
end

T1sel = T1r(:,selFlEx); T2sel = T2r(:,selFlEx);
x1 = mean(T1sel.Variables, 2); x1v = std(T1sel.Variables, [], 2);
x2 = mean(T2sel.Variables, 2); x2v = std(T2sel.Variables, [], 2);
x1t = filtfilt(filttrem,1,x1); x1g = filtfilt(filtgrasp,1,x1);
x2t = filtfilt(filttrem,1,x2); x2g = filtfilt(filtgrasp,1,x2);
x1te = envelope(x1t); x1ge = envelope(x1g);
x2te = envelope(x2t); x2ge = envelope(x2g);

T1wr = T1r(:,(end-2):end); T2wr = T2r(:,(end-2):end);
Y1 = T1wr.Variables; Y2 = T2wr.Variables;
Y1t = filtfilt(filttrem,1,Y1); Y2t = filtfilt(filttrem,1,Y2);
Y1te = envelope(Y1t); Y2te = envelope(Y2t);

figure; 
ax(1) = subplot(3,1,1);
patch([T1r.Time; flipud(T1r.Time)], [x1+x1v; flipud(x1-x1v)], 'b', 'FaceAlpha',.6, 'EdgeColor','none');
grid on; hold on; 
patch([T2r.Time; flipud(T2r.Time)], [x2+x2v; flipud(x2-x2v)], 'r', 'FaceAlpha',.4, 'EdgeColor','none');
plot(T1r.Time, x1, 'b', 'LineWidth',1.1); plot(T2r.Time, x2, 'r', 'LineWidth',1);
legend('1SD', '1SD', 'baseline', 'stim'); ylabel('Avg Flexion/Extension');
title('Grasping');
ax(2) = subplot(3,1,2);
plot(T1r.Time, x1t); grid on; hold on; plot(T2r.Time, x2t);
plot(T1r.Time, x1g); grid on; hold on; plot(T2r.Time, x2g);
legend('baseline tremor', 'stim tremor', 'baseline motor', 'stim motor'); 
ylabel('Avg Flexion/Extension');
ax(3) = subplot(3,1,3); 
plot(T1r.Time, Y1t); grid on; hold on; plot(T2r.Time, Y2t);
legend([T1wr.Properties.VariableNames, T2wr.Properties.VariableNames]);
ylabel('wrist position'); 
linkaxes(ax,'x');

[mean(x1te(T1i)), mean(x1ge(T1i))]
[mean(x2te(T2i)), mean(x2ge(T2i))]
[~,p1] = ttest2(x1te(T1i), x2te(T2i), "Vartype","unequal", "Tail","both")
[~,p2] = ttest2(x1ge(T1i), x2ge(T2i), "Vartype","unequal", "Tail","both")
mean(Y1te(T1i,:))
mean(Y2te(T2i,:))
[~,p3] = ttest2(Y1te(T1i,1), Y2te(T2i,1), "Vartype","unequal", "Tail","both")
[~,p4] = ttest2(Y1te(T1i,2), Y2te(T2i,2), "Vartype","unequal", "Tail","both")
[~,p5] = ttest2(Y1te(T1i,3), Y2te(T2i,3), "Vartype","unequal", "Tail","both")