% comparing constant vs dynamic AR parameters 

% parameters 
learnrate = 0.05;
ARord = 50; % # coefficients 
donorm = false;

% time points (sample) 
t1 = 10500; % start AR train 
t2 = t1+1000; % end AR train 
t3 = t2+ARord+1; % start eval AR (1)
t4 = t2+577; % start eval AR (2)
T = 100; % AR eval dur

%% file selection 
%[fn,fp] = uigetfile('*SegmentData*.mat', 'Choose Segmented Data File');
%SegmentedDataFullfile = fullfile(fp,fn);
SegmentedDataFullfile = '/Users/torenarginteanu/Documents/MATLAB/NeuralSystemID/AdaptAR/PD22N009_SegmentData_7e8979e58123db.mat';
load(SegmentedDataFullfile);
[fp,fn,fe] = fileparts(SegmentedDataFullfile);
subjname = upper(fn(1:8));
isEMU = contains(subjname, 'PY');

if isEMU
    freqrng = [4, 9]; % Hz (Theta band)
else
    freqrng = [13, 30]; % Hz (beta band)
end

% data selection 
tblBaseline = tblsBaseline{1};
chName = tblBaseline.Properties.VariableNames; 
chDesc = tblBaseline.Properties.VariableDescriptions;
chNameOnly = regexprep(upper(chName), '-REF', '');
chNameOnly = regexprep(chNameOnly, '\d+', ''); % no number
if isEMU
% only include macro channels both named (i.e. targeted) and identified for
% hippocampus with no uncertainty
chsel = ...
    (strcmp(chNameOnly,'LH')  | strcmp(chNameOnly,'RH')  | ... R/L hippo
     strcmp(chNameOnly,'LAH') | strcmp(chNameOnly,'RAH') | ... R/L ant hippo
     strcmp(chNameOnly,'LPH') | strcmp(chNameOnly,'RPH')) & ... R/L post hippo
    contains(lower(chDesc),'hip') & ~contains(chDesc,'?'); ...
%    & ~contains(chName,'BF'); 
else
    chsel = contains(chDesc, 'PrG'); % Precentral Gyrus
end
dtaBL = tblBaseline(:,chsel);
chselName = dtaBL.Properties.VariableNames;

% sample rate 
Fs = dtaBL.Properties.SampleRate;
if isnan(Fs)
    dt = seconds(diff(dtaBL.Time));
    if (min(dt) < 0) 
        error('Time series must be ascending.')
    end
    dtmean = mean(dt);
    if max(abs(dt-dtmean)) > .005
        error('Sample rate is not consistent.')
    end
    Fs = 1/dtmean;
end

%% subselect channels 

chsel1 = listdlg("PromptString","select channel", "SelectionMode","single", ...
    "ListString",chselName, ...
    "InitialValue", find(strcmp(chselName,channelNameRec)));
x = dtaBL{:,chsel1};
t = ((1:length(x))-1)/Fs; %t = t-t(t1);

chselname = dtaBL.Properties.VariableNames{chsel1}
chseldesc = dtaBL.Properties.VariableDescriptions{chsel1}

% filter 
BPF = fir1(1023, freqrng/(Fs/2));
x = filtfilt(BPF,1,x);

%% compute AR 

% constant 
x1 = x(t1:t2); 
mdl = ar(iddata(x1,[],1/Fs), ARord, 'yw');
x3c = myFastForecastAR(mdl, x(t2:t3), T);
x4c = myFastForecastAR(mdl, x(t2:t4), T);

% adaptive 
w1 = mdl.A; w1 = fliplr(-w1(2:end)/w1(1));
w = w1;
W = w;
for ti = t2:t3
    xi = x(ti); xxi = x((ti-ARord):(ti-1));
    w_ = updateWts(w, xxi, xi, learnrate, donorm);
    r = roots([1, -fliplr(w_)]);
    if max(abs(r)) < 1 % ensure stability
        w = w_;
    end
    W = [W; w];
end
w3 = w; mdl3 = [norm(w3)/norm(w1), -fliplr(w3)];
x3a = myFastForecastAR(mdl3, x(t2:t3), T);
for ti = (t3+1):t4
    xi = x(ti); xxi = x((ti-ARord):(ti-1));
    w_ = updateWts(w, xxi, xi, learnrate, donorm);
    r = roots([1, -fliplr(w_)]);
    if max(abs(r)) < 1 % ensure stability
        w = w_;
    end
    W = [W; w];
end
w4 = w; mdl4 = [norm(w4)/norm(w1), -fliplr(w4)];
x4a = myFastForecastAR(mdl4, x(t2:t4), T);

figure; 
subplot(2,1,1);
imagesc(diff(W,2)'); %colorbar; 
title('AR coeff change');
ylabel('coef index'); xlabel('sample');
subplot(2,1,2);
plot(rms(diff(W,2),2));
ylabel('mag diff');

%% plot time

clr = {[0.0660    0.4430    0.7450], ... blue 
       [0.8660    0.3290         0], ... red 
       [0.2310    0.6660    0.1960], ... green
       [0.5210    0.0860    0.8190], ... purple
       ...[0.6193    0.4627    0.0833], ... gold
       [0.6110    0.4660    0.1250], ... gold
       ...[0.9290    0.6940    0.1250], ... yellow
       ...[0.2588    0.4471    0.5294], ... teal
       [     0    0.6390    0.6390], ... teal
       [0.8190    0.0150    0.5450], ... pink
       [0.3720    0.1050    0.0310], ... brown
       [0.7170    0.1920    0.1720], ... dark red
       [0.0070    0.3450    0.0540], ... dark green
       [0.0620    0.2580    0.5010] ... dark blue 
       };
FaceAlpha = 0.6;

figure('Position',[1 1 1200 375], 'WindowStyle','normal', ...
    'Theme','light', 'Color','w'); 
plot(t, x, 'k', 'LineWidth',2); hold on; % grid on; 
plot(t(t3+(1:T)), x3c,      'color',clr{1}, 'LineWidth',3);
plot(t(t3+(1:T)), x3a, ':', 'color',clr{2}, 'LineWidth',3);
plot(t(t4+(1:T)), x4c,      'color',clr{1}, 'LineWidth',3);
plot(t(t4+(1:T)), x4a, ':', 'color',clr{2}, 'LineWidth',3);

xlim(t([t2-.5*T, t4+1.5*T]));
ax5 = gca(); ax5.FontSize = 14;
xlabel(' time', 'FontSize',16, ...
    'Units','Normalized', 'Position',[1,0.5], ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle'); 
xticks(t([t3 t4])); xticklabels({'t_1', 't_2'});
ylabel('Signal (\muV)', 'FontSize',16);
title('Model-Forecast Signal Example', 'FontSize',16);
subtitle(['Subject ',subjname], 'FontSize',16);
ax5.Box = false;
ax5.XAxisLocation = 'origin';
ax5.XAxis.TickLength = [.1 .025];
%ax5.XTick = unique(round(ax5.XTick)); % round to nearest second
%ax5.XTickLabels = {'1','2'};
ax5.XAxis.TickDirection = 'both';
ax5.YTick = [-50 0 50];
ax5.YAxis.TickDirection = 'both';

% time bar
yl = ylim; xl = xlim;
xbar = xl(1) + 0.01*diff(xl); xbar = xbar + [0 .25];
errorbar(mean(xbar), yl(2), 0,0, diff(xbar)/2, diff(xbar)/2, ...
    'k', 'LineWidth',2, 'CapSize',10);
text(mean(xbar), yl(2), '0.25 s', 'FontSize',14, ...
    'HorizontalAlignment','center', 'VerticalAlignment','top');

legend('Actual', 'Constant Model', 'Adaptive Model', ...
    'FontSize',16, 'location','southoutside', 'Orientation','horizontal');

%% plot coeffs before-after 

w1 = W(1,:); w2 = W(end,:); % 't1' and 't2' on the graph 
warrow = 10*(w2-w1);

% draw vertical arrowheads explicitly 
warrowpos = warrow; warrowneg = warrow;
warrowpos(warrowpos<0) = nan; 
warrowneg(warrowneg>0) = nan;
warrowpos = warrowpos + w1; warrowneg = warrowneg + w1;

% prepare figure 
figure('Position',[1 1 800 375], 'WindowStyle','normal', ...
    'Theme','light', 'Color','w'); 
tiledlayout('vertical', 'TileSpacing','tight'); % or none

% make small top accessory plot and bottom main plot
ax(1) = nexttile; ax(2) = nexttile([2,1]);

% graph coeff vals with direction of change
for curax = ax
axis(curax, 'tight');
stem(curax, w1, 'k', 'LineWidth',2); 
hold(curax, 'on'); quiver(curax, 1:length(w1),w1, zeros(size(w2)), warrow, 'off', ...
    'r', 'LineWidth',2, 'ShowArrowHead','off');
plot(curax, warrowpos, '^r', 'LineWidth',2, 'MarkerSize',2, 'MarkerFaceColor','auto'); 
plot(curax, warrowneg, 'vr', 'LineWidth',2, 'MarkerSize',2, 'MarkerFaceColor','auto');
curax.FontSize = 14;
end
linkaxes(ax, 'x');

% label/scale bottom main plot
axes(ax(2));
xlbl = xlabel('Tap #', 'FontSize',16);
ylabel('Coeff. Value', 'FontSize',16);
legend('t_1 value', 't_2 change (×10 scale)', 'Location','best', ...
    'FontSize',14);
yl = ylim; 
ybrk = max(w1(1:end-1)+warrow(1:end-1)); 
ybrk = ybrk + min(w1(1:end-1)+warrow(1:end-1)) - yl(1); 
ylim([yl(1), ybrk]);
xticks(length(w1))
xlblpos = xlbl.Position; ax(2).XAxisLocation = 'origin'; 
xlbl.Position = xlblpos; xlbl.Position(2) = yl(1);
ax(2).Box = false;
ax(2).XAxis.TickDirection = 'in';
ax(2).XRuler.TickLabelGapOffset = -30;

% label/scale top accessory plot
axes(ax(1));
title('Model Coefficient Update', 'FontSize',16)
ysz1 = ax(1).Position(4); ysz2 = ax(2).Position(4);
yscl = diff(ax(2).YLim)*ysz1/ysz2;
ycent = w1(end)+warrow(end)/2;
ylim(ycent + [-.5,.5]*yscl); 
yticks(round(w1(end), 2));
ax(1).Box = false;
ax(1).XAxis.Visible = false;

%% plot sys poles before/after 

% compute roots 
p1 = roots([1, -fliplr(w1)]); p2 = roots([1, -fliplr(w2)]);
parrow = (p2-p1); 

% configure inset showing greatest change
parrowlen = diag(parrow*parrow');
psel = parrowlen > 1e-4;
psel = [p1(psel); p2(psel)];
bselbox = [min(real(psel)), max(real(psel)); 
           min(imag(psel)), max(imag(psel))];
bselbox = bselbox + .25*[-1,1].*diff(bselbox,[],2);

% prepare figure
figure('Position',[100 100 875 375], 'WindowStyle','normal', ...
    'Theme','light', 'Color','w'); 
tiledlayout('horizontal', 'TileSpacing','compact'); 
ax2(1) = nexttile; ax2(2) = nexttile;

% make plot
for curax = ax2
axes(curax); axis(curax, 'tight');
viscircles(curax, [0,0],1, 'Color',.5*[1,1,1], 'LineStyle','--', 'LineWidth',3);
hold on; grid on;
plot(p1, 'ok', 'LineWidth',2); 
%plot(p2, 'or', 'LineWidth',1.5);
quiver(real(p1),imag(p1), real(parrow),imag(parrow), "off", ...
    'r', 'LineWidth',2, 'ShowArrowHead','on', 'MaxHeadSize',.5);
curax.Box = false;
curax.XAxisLocation = "origin";
curax.YAxisLocation = "origin";
rectangle(curax, 'Position',[bselbox(:,1); diff(bselbox,[],2)]', ...
    'EdgeColor','b', 'LineWidth',2)
curax.FontSize = 14;
xlabel('Re.', 'FontSize',14, 'FontWeight','bold', ...
    'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'Position',[0.1,0,-1]);
ylabel('Im.', 'FontSize',14, 'FontWeight','bold', ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
    'Position',[0,0.1,-1], 'Rotation',90);
end

% make inset plot
axes(ax2(2));
xlim(bselbox(1,:)); ylim(bselbox(2,:));
ax2(2).YAxisLocation = "right";
title(ax2(2),'Dynamic Poles', 'FontSize',16)

% annotate main plot
title(ax2(1),'All Model Poles', 'FontSize',16)
lgd = legend(ax2(1), 't_1 Value', 't_2 Change', ...
    'FontSize',14, 'location','westoutside');
lgd.Layout.Tile = 'west';

%% helper(s)

    % Update based on gradient wrt weights; no guarantee of stability.
    function w = updateWts(w, x, y, stepsize, donorm)
        ypred = w*x;
        E = y-ypred; del = x*E;
        if donorm
            del = del./(x'*x + eps);
        end
        w = w + stepsize*del';
    end