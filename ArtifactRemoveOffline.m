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
tblsName = {'tblsTrig'; 'tblsStimNoTrig'; 'tblsSrl'; 'tblsMisc'};
tblsDesc = tblsName;
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
tblsDesc = tblsDesc(selidx); tblsList = tblsList(selidx); tblsName = tblsName(selidx);
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

tblsListOut = cell(size(tblsList));

for Li = 1:height(tblsList)
    tbls = tblsList{Li};
    tblsOut = cell(height(tbls),1);
    for Ti = 1:height(tbls)
        dta = tbls{Ti,1}; % only process the "main table"
        dta = dta(:, chanlistsel);
        dta = sortrows(dta, 'Time'); % order time ascending, if not already done
        dta = dta - DCOSBL; % correct DC offset 
        t = dta.Time; ts = seconds(dta.Properties.TimeStep);
        if isnan(ts)
            ts = mode(seconds(diff(t)));
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
noises = nan(2*N, width(dta), length(noise_inds), 2);
for ind = 1:size(noises, 3)
    dta_inds = dta{noise_inds{ind},:}; dtaKal_inds = dtaKal{noise_inds{ind},:};
    dta_inds = dta_inds(1:min(height(dta_inds), height(noises)), :);
    dtaKal_inds = dtaKal_inds(1:height(dta_inds),:);
    noises(1:height(dta_inds),:,ind,1) = dta_inds;
    noises(1:height(dta_inds),:,ind,2) = dtaKal_inds;
end
clear noise_ind ind ind_curr noise_ind_new dta_inds
artavg = mean(noises,3,'omitnan'); % avg waveform accross stims 
artstd = std(noises,[],3,'omitnan'); % var in waveform across trials
artdif = noises - artavg; artdifrsq = sqrt(sum(artdif.^2,1,'omitnan'));
bestartdif = diff(artdifrsq,[],4); bestartdif = sum(bestartdif.^2,3,'omitnan');
bestartdif = squeeze(bestartdif); [~,bestartdif] = sort(bestartdif);
artdifmin = min(artdif,[],1); % most neg val each trial 
artdifmax = max(artdif,[],1); % most pos val each trial
artdifavg = mean(artdif,1,'omitnan'); 
artimbal = (artdifrsq - artdifavg)./(artdifrsq+eps); % polarity imbalance each trial

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
    chanbestprms(1), chanbestprms(end), ...
    bestSNRupper(1), bestSNRupper(end), ...
    bestSNRlower(1), bestSNRlower(end), ...
    bestartdif(1), bestartdif(end)]; 
chandisp = unique(chandisp); H = length(chandisp);
chandispname = chanlistsel(chandisp);
if ~sum(strcmp(channelNameRec, chandispname))
    if sum(strcmp(channelNameRec, chanlistsel))
        chandispname = [chandispname, channelNameRec];
    end
end
for chnmst = channelNameStim
    if ~sum(strcmp(chnmst, chandispname))
        if sum(strcmp(chnmst, chanlistsel))
            chandispname = [chandispname, chnmst];
        end
    end
end
figure('Units','normalized', 'Position',[.05,.05,.9,.9]); 
for ch = 1:length(chandispname)
    ax(ch) = subplot(H+1,1,ch);
    plot(dta, chandispname(ch)); grid on; hold on; 
    plot(dtaKal, chandispname(ch));
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
for ind = 1:size(noises,3)
    yi = noises(:,chandisp(ch),ind,1); % with artifact
    xi = noises(:,chandisp(ch),ind,1); % artifact-removed
    ti = 1:length(yi); ti = (ti-1)*ts;
    plot(ti, yi, ...
        'Color',colorwheel(ind/length(noises)), ...
        'LineWidth',1.5-ind/length(noises));
end
ylabel(chandispname(ch)); xlabel('time from stim (s)');
title('artifacts')
end
linkaxes(ax2, 'x');

pause(.001); drawnow; pause(.001);

% show characteristic artifact waveform and removal
NumTraceToShow = 32; % max
tracetypes = cat(5, artdifrsq, artdifmin, artdifmax, artimbal);
NumTraceToShow = floor(NumTraceToShow/(2*size(tracetypes,5)));
tracetypes = tracetypes(:,chandisp,:,1,:); % selected channels, with artifact
tracetypes = sum(tracetypes.^2,2,'omitnan'); % sum all selected channels
tracetypes = squeeze(tracetypes); % trial x type 
tracebest = zeros(size(tracetypes));
for itype = 1:size(tracetypes,2)
    [~,tracebest(:,itype)] = sort(tracetypes(:,itype));
end
showtraceinds = [tracebest(1:NumTraceToShow,:); tracebest((end-NumTraceToShow+1):end,:)];
showtraceinds = showtraceinds(:); showtraceinds = unique(showtraceinds);
figure; 
for ch = 1:H
ax3(ch) = subplot(H,1,ch);
yavg = artavg(:,chandisp(ch),1,1); xavg = artavg(:,chandisp(ch),1,2);
ystd = artstd(:,chandisp(ch),1,1); xstd = artstd(:,chandisp(ch),1,2);
patch([ti,fliplr(ti)], [yavg+ystd;flipud(yavg-ystd)], 'b', 'FaceAlpha',.5, 'EdgeColor','none');
hold on; grid on;
patch([ti,fliplr(ti)], [xavg+xstd;flipud(xavg-xstd)], 'r', 'FaceAlpha',.5, 'EdgeColor','none');
for ind = 1:length(showtraceinds)
    yi = noises(:,chandisp(ch),showtraceinds(ind),1);
    xi = noises(:,chandisp(ch),showtraceinds(ind),2);
    plot(ti,yi, ':b', "LineWidth",.75);
    plot(ti,xi, ':r', "LineWidth",.75);
end
plot(ti,yavg, "Color",'b', "LineWidth",2);
plot(ti,xavg, "Color",'r', "LineWidth",1.7);
ylabel(chandispname(ch)); xlabel('time from stim (s)');
title('artifacts')
end
linkaxes(ax3, 'x');

pause(.001); drawnow; pause(.001);

dtaKal = dtaKal + DCOSBL; % replace DC offset 
tblsOut{Ti,1} = dtaKal;

    end
    tblsListOut{Li} = tblsOut;
end

%% saving 

% assign new variables
for Li = 1:height(tblsListOut)
    tblsName{Li} = [tblsName{Li},'_ArtifactRemoved'];
    eval([tblsName{Li},' = tblsListOut{Li};']);
end

% stamp saved output with version
thisfilever = getFileVersion(thisfilename);
[~,thisfilename] = fileparts(thisfilename);
svname = [fn,'_',thisfilename,'_',thisfilever];

% save vars (table lists)
save(fullfile(fp,svname), tblsName, 'chanlistsel', "-v7.3");