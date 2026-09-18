% comparing constant vs dynamic AR parameters 

% parameters 
learnrate = 0.05;
ARord = 50; % # coefficients 
donorm = false;

% time points (sample) 
t1 = 1; % start AR train 
t2 = 1000; % end AR train 
t3 = 1060; % start eval AR (1)
t4 = 1500; % start eval AR (2)
T = 250; % AR eval dur

%% file selection 
[fn,fp] = uigetfile('*SegmentData*.mat', 'Choose Segmented Data File');
SegmentedDataFullfile = fullfile(fp,fn);
load(SegmentedDataFullfile);
[~,fn,fe] = fileparts(fn);
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
    "ListString",chselName);
x = dtaBL{:,chsel1};
t = ((1:length(x))-1)/Fs;

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
for ti = t2:t3
    xi = x(ti); xxi = x((ti-ARord):(ti-1));
    w_ = updateWts(w, xxi, xi, learnrate, donorm);
    r = roots([1, -fliplr(w_)]);
    if max(abs(r)) < 1 % ensure stability
        w = w_;
    end
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
end
w4 = w; mdl4 = [norm(w4)/norm(w1), -fliplr(w4)];
x4a = myFastForecastAR(mdl4, x(t2:t4), T);

%% plot

figure; 
plot(t, x, 'w', 'LineWidth',1.5); hold on; grid on; 
plot(t(t3+(1:T)), x3c, 'b', 'LineWidth',1);
plot(t(t3+(1:T)), x3a, '--r', 'LineWidth',1);
plot(t(t4+(1:T)), x4c, 'b', 'LineWidth',1);
plot(t(t4+(1:T)), x4a, '--r', 'LineWidth',1);
xlim(t([t2-.5*T, t4+1.5*T]));

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