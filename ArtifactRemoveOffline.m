%% Offline Artifact Removal 
% This code is designed to remove all stimulus artifacts from data that has
% been segmented so that data with stimulation can be further analyzed.
% This should use the absolute best-performing offline method of artifact
% removal, regardless of causality or real-time applicability. This is not
% intended to compare or evaluate different artifact removal techniques. 

%% obtain segmented data 

% file selection
thisfilename = mfilename("fullpath");
[fn,fp] = uigetfile('*SegmentData*.mat', 'Choose Segmented Data File');
SegmentedDataFullfile = fullfile(fp,fn);
load(SegmentedDataFullfile);
[~,fn,fe] = fileparts(fn);

% display details 
disp("Stimulus on channel(s): "); disp(channelNameStim);
disp("Trigger on channel(s): "); disp(channelNameTrig);

%% user selects which data to artifact-remove

% select table lists of interest 
tblsDesc = {'tblsTrig'; 'tblsStimNoTrig'; 'tblsSrl'; 'tblsMisc'};
tblsList = cellfun(@eval, tblsDesc, 'UniformOutput', false);
for Li = 1:height(tblsDesc)
    TC = tblsList{Li}; TCdescs = '';
    for Ti = 1:height(TC)
        T = TC{Ti};
        Tdesc = T.Properties.Description; Tdesc = char(Tdesc);
        if ~strcmpi(TCdescs, Tdesc)
            TCdescs = [TCdescs,' & ',Tdesc];
        end
    end
    TCdescs = TCdescs(4:end); % remove initial ' & '
    tblsDesc{Li} = [tblsDesc{Li},': ',TCdescs];
end
[selidx,selmade] = listdlg("ListString",tblsDesc, ...
    "SelectionMode","multiple", "PromptString","Select Data With Artifact", ...
    "ListSize",[500,300]);
if ~selmade
    error('Selection must be made.')
end
tblsDesc = tblsDesc(selidx); tblsList = tblsList(selidx);
clear Li Ti T TC TCdescs Tdesc selidx selmade

% subselect tables of interest 
for Li = 1:height(tblsList)
    TC = tblsList{Li}; 
    if height(TC) > 1
        TC1 = TC(:,1);
        TCdescs = cellfun(@(T) T.Properties.Description, TC1, 'UniformOutput',false);
        [selidx,selmade] = listdlg("ListString",TCdescs, ...
            "SelectionMode","multiple", "PromptString","Select Data With Artifact", ...
            "ListSize",[500,300]);
        if ~selmade
            error('Selection must be made.')
        end
        TC = TC(selidx,:);
        tblsList{Li} = TC;
    end
end
clear Li TC TC1 TCdescs selidx selmade

% select channels of interest 
dtaBL = tblsBaseline{1,1};
chanlist = dtaBL.Properties.VariableNames; chanlist = string(chanlist);
chandesc = dtaBL.Properties.VariableDescriptions; chandesc = string(chandesc);
[chanselidx,chanselmade] = listdlg("ListString", chanlist+": "+chandesc, ...
    "SelectionMode","multiple", "PromptString","Select channels to INCLUDE", ...
    "ListSize",[500,300]);
chanlistsel = chanlist(chanselidx);

%% train baseline model 

if chanselmade
    dtaBL = dtaBL(:, chanlistsel);
end

% remove DC offset 
DCOSBL = mean(dtaBL);
dtaBL = dtaBL - DCOSBL;

[predBL, ~, ~, tstEval, A] = fitLTIauton(dtaBL);
tstEval

% use training residual to estimate process noise Q
res = predBL - dtaBL; % residual 
Qcov = cov(res.Variables);
chanrms = rms(res.Variables); chanprms = chanrms./rms(dtaBL.Variables);
[~,chanbestprms] = sort(chanprms);

%% loop through table list 

for Li = 1:height(tblsList)
    tbls = tblsList{Li};
    for Ti = 1:height(tbls)
        dta = tbls{Ti,1}; % only process the "main table"
        dta = dta(:, chanlistsel);
        dta = sortrows(dta, 'Time'); % order time ascending, if not already done
        dta = dta - DCOSBL; % correct DC offset 
        t = dta.Time; ts = seconds(dta.Properties.TimeStep);
        if isnan(ts)
            ts = mode(seconds(diff(ts)));
        end

        % determine trigger (noise ref) signal, g
        evtbl = dta.Properties.Events;
        evlbl = unique(evtbl.EventLabels);
        [selidx,selmade] = listdlg("ListString",evlbl, ...
            "SelectionMode","single", "PromptString","Select Noise Reference", ...
            "ListSize",[500,300]);
        if ~selmade
            error('Selection must be made.')
        end
        gname = evlbl(selidx);
        selrow = strcmp(evtbl.EventLabels, gname);
        gtbl = evtbl(selrow,:);
        gt = gtbl.Time; 
        g = zeros(size(t));
        for gti = gt'
            [gtidiff,gi] = min(abs(seconds(t-gti)));
            if gtidiff < 5*ts % tolerance; otherwise ignored
                g(gi) = 1;
            end
        end
        clear gi gti gtidiff selmade selidx selrow

% LMS setup and pretraining 
N = 50; % filter # taps
stepsize = .5; % learn rate for gradient descent 
W0 = preTrainWtsLMS(g,dta,N,4,false); % "optimal" LMS weights 

% Kalman filter with no adaptive state estimate 
dtaKal = AdaptKalmanAuton(g,dta,[],0,A,stepsize,Qcov,N,W0,true,false,true);

% characterize artifact waveform 
noise_ind = find(g); 
noise_ind_new = diff(noise_ind) > 2;
noise_ind_new = find(noise_ind_new) + 1;
noise_ind_new = noise_ind(noise_ind_new);
noise_inds = {};
for ind = 2:length(noise_ind_new)
    ind_curr = noise_ind_new(ind-1) : (noise_ind_new(ind)-1);
    noise_inds = [noise_inds, (ind_curr)];
end
clear noise_ind ind ind_curr noise_ind_new
noise_tbls = cellfun(@(inds) dta(inds,:), noise_inds, 'UniformOutput',false);

% assess noise level by channel 
SNRupper = mag2db(rms(dta.Variables)./(rms(dta.Variables)-rms(dtaBL.Variables)));
[~,bestSNRupper] = sort(SNRupper);
SNRlower = nan(size(SNRupper));
for ch = 1:width(dta)
    SNRlower(ch) = snr(dta{:,ch}, dta{:,ch}-dtaKal{:,ch});
end
[~,bestSNRlower] = sort(SNRlower);

% display some example channels 
chandisp = [...
    chanbestprms(1); chanbestprms(end); ...
    bestSNRupper(1); bestSNRupper(end); ...
    bestSNRlower(1); bestSNRlower(end)]; 
chandisp = unique(chandisp); H = height(chandisp);
chandisp = chanlistsel(chandisp);
if ~sum(strcmp(channelNameRec, chandisp))
    if sum(strcmp(channelNameRec, chanlistsel))
        chandisp = [chandisp, channelNameRec];
    end
end
for chnmst = channelNameStim
    if ~sum(strcmp(chnmst, chandisp))
        if sum(strcmp(chnmst, chanlistsel))
            chandisp = [chandisp, chnmst];
        end
    end
end
figure('Units','normalized', 'Position',[.05,.05,.9,.9]); 
for ch = 1:length(chandisp)
    ax(ch) = subplot(H+1,1,ch);
    plot(dta, chandisp(ch)); grid on; hold on; 
    plot(dtaKal, chandisp(ch));
end
ax(ch+1) = subplot(H+1,1,ch+1);
plot(t, g); grid on;
linkaxes(ax, 'x');

pause(.001); drawnow; pause(.001);

% show overlapping artifact
figure; 
for ch = 1:H
ax2(ch) = subplot(H,1,ch);
hold on; grid on;
for TBLi = 1:length(noise_tbls)
    tbl = noise_tbls{TBLi};
    tbl.Time = tbl.Time - tbl.Time(1);
    plot(seconds(tbl.Time), tbl.(chandisp(ch)), ...
        'Color',colorwheel(TBLi/length(noise_tbls)), ...
        'LineWidth',1.5-TBLi/length(noise_tbls));
end
ylabel(chandisp(ch)); xlabel('time (s)');
title('artifacts')
end
linkaxes(ax2, 'x');
xlim([0 2*N*ts])

pause(.001); drawnow; pause(.001);

    end
end