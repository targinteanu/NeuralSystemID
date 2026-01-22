%% Feature Selection Tool
% Explore preprocessing steps to see which features depend on (A) previous
% values of the same features and (B) stimulation 

hzns = [.05]; % s

thisfilename = mfilename("fullpath");

% init empty tbl lists 
tblsTrig_ArtifactRemoved = {};
tblsTrig = {};
tblsStimNoTrig_ArtifactRemoved = {};
tblsStimNoTrig = {};

ViewSpectrum

%% organize into tables 

featnames = ["Raw", "Filtered", "Hilbert", "MagPhase"];

alltbls = cell(4,length(featnames));
% stim type {baseline, other non-stim, cort stim, depth stim} x feature
% all wrapped in cells so there can be multiple 
alltbls{1,1} = {dtaBL};

% add "misc" tables
tblsToOrganize = tblsMisc(:,1); 
if chanselmade
    for Tj = 1:height(tblsToOrganize)
        tblsToOrganize{Tj} = tblsToOrganize{Tj}(:,chanselidx);
    end
end

% collect by stim type 
tblsToOrganize = [...
    tblsToOrganize; ... misc, selected channels only
    selectTbls(tblsTrig, tblsTrig_ArtifactRemoved, chanselmade, chanlistsel, chanselidx); ...
    selectTbls(tblsStimNoTrig, tblsStimNoTrig_ArtifactRemoved, chanselmade, chanlistsel, chanselidx)];
tblsToOrganizeDescs = cellfun(@(T) string(T.Properties.Description), tblsToOrganize);

% other non-stim 
[selidx,selmade] = listdlg("ListString",tblsToOrganizeDescs, ...
    "SelectionMode","multiple", "ListSize",[500,300], ...
    "PromptString","Select NON-stimulation");
if selmade
    alltbls{2,1} = tblsToOrganize(selidx);
end

% cort stim 
[selidx,selmade] = listdlg("ListString",tblsToOrganizeDescs, ...
    "SelectionMode","multiple", "ListSize",[500,300], ...
    "PromptString","Select CORTICAL stimulation");
if selmade
    alltbls{3,1} = tblsToOrganize(selidx);
end

% depth stim 
[selidx,selmade] = listdlg("ListString",tblsToOrganizeDescs, ...
    "SelectionMode","multiple", "ListSize",[500,300], ...
    "PromptString","Select DEPTH stimulation");
if selmade
    alltbls{4,1} = tblsToOrganize(selidx);
end

% TO DO: if any of the above tables has missing or nan for a time period
% greater than tol, break it into cell array of tables 

% homogenize variable descriptions and units 
for H = 1:height(alltbls)
    tbls = alltbls{H,1};
    for h = 1:height(tbls)
        tbl = tbls{h};
        if isempty(tbl.Properties.VariableUnits)
            tbl.Properties.VariableUnits = repmat("", 1, width(tbl));
        end
        if isempty(tbl.Properties.VariableDescriptions)
            tbl.Properties.VariableDescriptions = repmat("", 1, width(tbl));
        end
        tbls{h} = tbl;
    end
    alltbls{H,1} = tbls;
end

%% processing steps for all conditions 

% notch out power line noise and harmonics
f0 = 60; % power line fundamental 
qFactor = 35;
notchB = 1; notchA = 1;
for h = f0:f0:(SampleRate/2)
    % add a notch at harmonic h
    [notchBh,notchAh] = iirnotch(h/(SampleRate/2), (h/(SampleRate/2))/qFactor);
    notchB = conv(notchB, notchBh); notchA = conv(notchA, notchAh);
end
for Ti = 1:height(alltbls)
    disp(['Notch Filtering: ',num2str(Ti),' of ',num2str(height(alltbls))])
    for Tj = 1:height(alltbls{Ti,1})
        T = alltbls{Ti,1}{Tj};
        for c = 1:width(T)
            T{:,c} = filtfilt(notchB,notchA,T{:,c});
        end
        alltbls{Ti,1}{Tj} = T;
    end
end

%%{
% reref to common average 
for Ti = 1:height(alltbls)
    disp(['Rereferencing: ',num2str(Ti),' of ',num2str(height(alltbls))])
    for Tj = 1:height(alltbls{Ti,1})
        T = alltbls{Ti,1}{Tj};
        T.Variables = T.Variables - mean(T.Variables, 2, 'omitnan');
        alltbls{Ti,1}{Tj} = T;
    end
end
%}

% visualize results 
sigBL = alltbls{1,1}{1};
pause(.5); fig_spec2 = figure; periodogram(sigBL.Variables, [], [], SampleRate, 'power');
pause(.001); drawnow; pause(.001);

%% calc features 

% filter 
for Ti = 1:height(alltbls)
    disp(['Band Filtering: ',num2str(Ti),' of ',num2str(height(alltbls))])
    TTraw = alltbls{Ti,1};
    TTfilt = cell(size(TTraw));
    for Tj = 1:height(TTraw)
        Traw = TTraw{Tj};
        Tfilt = [];

        %{
        % broad band
            bpf = buildFIRBPF(SampleRate, bandbounds(1), bandbounds(end));
            varnames = string(Traw.Properties.VariableNames)+" all";
            varunits = string(Traw.Properties.VariableUnits);
            vardescs = string(Traw.Properties.VariableDescriptions);
            Xfiltall = filtfilt(bpf,1,Traw.Variables);
            Tfiltall = array2timetable(Xfiltall,...
                "RowTimes",Traw.Time, "VariableNames",varnames);
                Tfiltall.Properties.VariableUnits = varunits;
                Tfiltall.Properties.VariableDescriptions = vardescs;
            Tfiltall.Properties.Events = Traw.Properties.Events;
            Xfiltall = envelope(Xfiltall);
        %}

        % specific bands 
        for b = 1:length(bandnames)
            bpf = buildFIRBPF(SampleRate, bandbounds(b), bandbounds(b+1), 2, 201);
            varnames = string(Traw.Properties.VariableNames)+" "+bandnames(b);
            vardescs = string(Traw.Properties.VariableDescriptions);
            varunits = string(Traw.Properties.VariableUnits);
            tic
            Xfilt = filtfilt(bpf,1,Traw.Variables);
            timerdur = toc;
            disp(" - filter band "+bandnames(b)+" completed after "+...
                string(num2str(timerdur))+" seconds.")
            %Xfilt = Xfilt ./ (Xfiltall + eps); % unitless
            Tfilt_ = array2timetable(Xfilt,...
                "RowTimes",Traw.Time, "VariableNames",varnames);
                Tfilt_.Properties.VariableDescriptions = vardescs;
                Tfilt_.Properties.VariableUnits = varunits;
            Tfilt = [Tfilt, Tfilt_];
        end

        %TTfilt{Tj} = [Tfilt, Tfiltall];
        Tfilt.Properties.Events = Traw.Properties.Events;
        TTfilt{Tj} = Tfilt;
        TTfilt{Tj}.Properties.Description = Traw.Properties.Description;
    end
    alltbls{Ti,2} = TTfilt;
end

%% hilbert, mag/phase

% pink noise correction factor 
%{
% pinkP = (10.^pinkcoef(2,:)) .* bandcent'.^pinkcoef(1,:); 
f = [log10(bandcent); ones(size(bandcent))]';
pinkP = f*pinkcoef;
% row = freq band; col = channel
pinkP = sqrt( 10.^(pinkP/10)' ); % Power (dB) -> amplitude
%pinkP = pinkP'; 
pinkP = pinkP(:)'; % unfolded 
%}
Tfilt = alltbls{1,2}{1};
Xabs = envelope(Tfilt.Variables);
P = log10(Xabs.^2)';
F = repmat(bandcent, length(chanlistsel), 1); F = F(:);
%F = repmat(F, 1, width(P));
%F = F(:); P = P(:);
F = [log10(F), ones(size(F))];
pinkcoef2 = F\P;
%pinkP = [log10(bandcent)', ones(size(bandcent))']*pinkcoef2;
pinkP = F*pinkcoef2;
pinkP = sqrt(10.^pinkP');
%pinkP = repmat(pinkP, length(chanlistsel), 1); pinkP = pinkP(:)';
pinkP = median(pinkP);

for Ti = 1:height(alltbls)
    disp(['Hilbert: ',num2str(Ti),' of ',num2str(height(alltbls))])
    TTfilt = alltbls{Ti,2};
    TThilb = cell(size(TTfilt)); TTmaph = cell(size(TTfilt));
    for Tj = 1:height(TTfilt)
        Tfilt = TTfilt{Tj};
        varnames = string(Tfilt.Properties.VariableNames);
        varunits = string(Tfilt.Properties.VariableUnits);
        vardescs = string(Tfilt.Properties.VariableDescriptions);
        Xhilb = hilbert(Tfilt.Variables);
        Xabs = abs(Xhilb);
        
        Xhilb = Xhilb./pinkP; Xabs = abs(Xhilb); % corrected
        %{
        TTpink = array2timetable(pinkcoef', "RowTimes",Tfilt.Time, ...
            "VariableNames",["Pink1", "Pink2"]); 
        TTpink.Properties.VariableUnits = {'', 'uV^2/s'};
        %}

        TThilb{Tj} = [...
            array2timetable([real(Xhilb)], "RowTimes",Tfilt.Time, ...
                "VariableNames",varnames+" real"), ...
            array2timetable([imag(Xhilb)], "RowTimes",Tfilt.Time, ...
                "VariableNames",varnames+" imag")];
        TThilb{Tj}.Properties.VariableUnits = [varunits, varunits];
        TThilb{Tj}.Properties.VariableDescriptions = [vardescs, vardescs];
        %TThilb{Tj} = [TThilb{Tj}, TTpink];
        TThilb{Tj}.Properties.Description = Tfilt.Properties.Description;
        TThilb{Tj}.Properties.Events = Tfilt.Properties.Events;

        TTmaph{Tj} = [...
            array2timetable([Xabs], "RowTimes",Tfilt.Time, ...
                "VariableNames",varnames+" mag"), ...
            array2timetable([angle(Xhilb)], "RowTimes",Tfilt.Time, ...
                "VariableNames",varnames+" phase")];
        TTmaph{Tj}.Properties.VariableUnits = [varunits, ...
            repmat("radians",1, width( Xhilb ))];
        TTmaph{Tj}.Properties.VariableDescriptions = [vardescs, vardescs];
        %TTmaph{Tj} = [TTmaph{Tj}, TTpink];
        TTmaph{Tj}.Properties.Description = Tfilt.Properties.Description;
        TTmaph{Tj}.Properties.Events = Tfilt.Properties.Events;

    end
    alltbls{Ti,3} = TThilb; alltbls{Ti,4} = TTmaph;
end

%% evaluate difference between stim and baseline 

pause(.5); fig_stimresponse = figure('Units','normalized', 'Position',[.05,.05,.5,.9]); 
pause(.01); drawnow; pause(.01);
for feat = 1:size(alltbls,2)
    tblsBL = alltbls{1,feat};
    T = [];
    for Tj = 1:height(tblsBL)
        T = [T; tblsBL{Tj}];
    end
    XBL = T.Variables;
    XBL = single(XBL); % single precision to save runtime and memory
    varnames = T.Properties.VariableNames;
    subplot(size(alltbls,2),1,feat);
    stimnames = repmat("", 1, size(alltbls,1));
    for stim = 2:size(alltbls,1)
        %subplot(size(alltbls,1)-1, size(alltbls,2), size(alltbls,2)*(stim-2)+feat);
        tbls = alltbls{stim,feat};
        T = [];
        for Tj = 1:height(tbls)
            T = [T; tbls{Tj}];
        end
        if numel(T)
        stimnames(stim) = string(T.Properties.Description);
        X = T.Variables;
        X = single(X); % single precision to save runtime and memory
        p = nan(1,width(X));
        for c = 1:width(X)
            [~,~,p(c)] = kstest2(XBL(:,c), X(:,c));
        end
        stem((p)); hold on; grid on; 
        end
    end
        if length(p) > length(chanselidx)
            xticks(1:length(chanselidx):length(p)); 
            xticklabels(varnames(1:length(chanselidx):length(p)));
        else
            xticks(1:length(p)); xticklabels(varnames);
        end
        %{
        if stim == 2
            title(featnames(feat));
        end
        if feat == 1
            ylabel({char(T.Properties.Description), 'KS vs baseline'});
        end
        if stim == size(alltbls,1)
            xlabel('channel/measurement');
        end
        %}
        xlabel('channel/measurement'); ylabel('KS vs baseline');
        legend(stimnames(2:end), 'Location','best')
        pause(.001); drawnow; pause(.001);
    %end
end

%% organize into regularly spaced arrays for each hzn
% TO DO: disp output progress 

hznsN = ceil(hzns*SampleRate); % samples 

allarr = cell([size(alltbls), length(hzns)]);
for feat = 1:size(allarr,2)
    for stim = 1:size(allarr,1)
        tbls = alltbls{stim,feat};
        for hzn = 1:size(allarr,3)
            disp(['Dataset: ',num2str((feat-1)*size(allarr,2)+stim),' of ',num2str(numel(alltbls)),...
                '; horizon ',num2str(hzn),' of ',num2str(size(allarr,3))])
            inp = []; oup = [];
            inp = single(inp); oup = single(oup); % single precision to save runtime and memory
            for Tj = 1:height(tbls)
                % TO DO: this assumes all tbls columns are same order!
                X = tbls{Tj}.Variables; % assume reg spaced time 
                X = single(X); % single precision to save runtime and memory
                %for startind = 0:(hznsN(hzn)-1)
                startind = 0;
                    inpoup = X((1:hznsN(hzn):end)+startind, :);
                    inp = [inp; inpoup(1:(end-1),:)];
                    oup = [oup; inpoup(2:end,:)];
                %end
            end
            allarr{stim,feat,hzn} = cat(3,inp,oup);
        end
    end
end

%% evaluate "autonomous" relationship of "states" 

nbin = 100;
fig_auton = zeros(1,size(allarr,3));

Amin = inf; Amax = -inf; imgs = cell(size(allarr));
for hzn = 1:size(allarr,3)
    spind = 1;
    pause(.5); fig_auton(hzn) = figure('Units','normalized', 'Position',[.05,.05,.9,.9]); 
    sgtitle([num2str(1000*hzns(hzn)),' ms prediction']);
    pause(.01); drawnow; pause(.01);
    for stim = 1:size(allarr,1)
        for feat = 1:size(allarr,2)
            if numel(allarr{stim,feat,hzn})
            varnames = alltbls{stim,feat}{1}.Properties.VariableNames;
            imgs{stim,feat,hzn} = subplot(size(allarr,1), size(allarr,2), spind); spind = spind+1;
            inp = allarr{stim,feat,hzn}(:,:,1);
            oup = allarr{stim,feat,hzn}(:,:,2);
            %A = (inp'*inp)^-1*inp'*oup;
            A = nan(width(inp),width(oup));
            for chinp = 1:height(A)
                for choup = 1:chinp
                    a = DistanceCorrelation(inp(:,chinp), oup(:,choup));
                    A(chinp,choup) = a; A(choup,chinp) = a;
                end
            end
            Amin = min(Amin, min(A(:))); Amax = max(Amax, max(A(:)));
            imagesc(A, [0 1]); colorbar;

            if width(A) > length(chanselidx)
                xticks(1:length(chanselidx):width(A)); 
                xticklabels(varnames(1:length(chanselidx):width(A)));
                yticks(1:length(chanselidx):height(A)); 
                yticklabels(varnames(1:length(chanselidx):height(A)));
            else
                xticks(1:width(A)); xticklabels(varnames);
                yticks(1:height(A)); yticklabels(varnames);
            end

            if stim == 1
                % label col with feat type
                title(featnames(feat));
            end
            if feat == 1
                % label row with stim type 
                %if stim == 1
                %    ylabel('Baseline');
                %else
                    ylabel(alltbls{stim,feat}{1}.Properties.Description)
                %end
            end

            pause(.001); drawnow; pause(.001);
            end
        end
    end
end
%{
for hzn = 1:size(allarr,3)
    for stim = 1:size(allarr,1)
        for feat = 1:size(allarr,2)
            clim(imgs{stim,feat,hzn}, [Amin, Amax]);
        end
    end
end
%}

%% simplify all tables for saving 

selfeat = listdlg("PromptString","Select feature(s) to save", ...
    "ListString",featnames, "SelectionMode","multiple");
selfeatname = featnames(selfeat);
svtbls = alltbls(:,selfeat); 

t0 = svtbls{1}{1}.Time(1); 

%sizethresh = 2^34; % max variable size (bytes) to save in single document

% stack each stim type into one table 
for W = 1:width(svtbls)
for H = 1:height(svtbls)
    TBL = [];
    tbls = svtbls{H,W};
    for h = 1:height(tbls)
        tbl = tbls{h};
        if numel(tbl)
            tbl.Time = tbl.Time - t0;
        end
        if numel(tbl.Properties.Events)
            tbl.Properties.Events = eventtable(...
                tbl.Properties.Events.Time - t0, ...
                "EventLabels",tbl.Properties.Events.EventLabels, ...
                "EventEnds",tbl.Properties.Events.EventEnds - t0);
        else
            tbl.Properties.Events = [];
        end
        if numel(TBL)
            %sz = whos('TBL'); sz = sz.bytes;
            TBL = [TBL; tbl];
        else
            TBL = tbl;
        end
    end
    svtbls{H,W} = TBL;
end
end

%% saving 
disp('Saving preprocessed data...')

thisfilever = getFileVersion(thisfilename);
[~,thisfilename] = fileparts(thisfilename);
svname = [thisfilename,'_',thisfilever];

% make save folder 
fp = [fp,svname];
mkdir(fp);
svname = [fn,'_',svname]; % tag with names/versions of all previous steps
svname = fullfile(fp, svname);

% save each table 
for W = 1:width(svtbls)
for H = 1:height(svtbls)
    TT = svtbls{H,W}; 

    % data table 
    if numel(TT)
    svname_ = string(svname)+" - "+selfeatname(W)+" - "+string(TT.Properties.Description);
    disp(" - "+svname_)

    EV = TT.Properties.Events;
    SvData = [(single(seconds(TT.Time))), (single(TT.Variables))];
    SvData = array2table(SvData, ...
        "VariableNames", ['Time (s)', TT.Properties.VariableNames]);
    writetable(SvData, svname_+" - data.csv");

    % info 
    if contains(selfeatname, "Hilbert") || contains(selfeatname, "MagPhase")
        SvIfoPink = arrayfun(@num2str, pinkP, 'UniformOutput', false);
        SvIfoPink = ['Pink Correction', SvIfoPink, SvIfoPink];
    else
        SvIfoPink = [];
    end
    SvIfo = ['Channel Name', TT.Properties.VariableNames; ...
             'Unit', TT.Properties.VariableUnits; ...
             'Description', TT.Properties.VariableDescriptions; ...
             SvIfoPink]';
    SvIfo{1, width(SvIfo)+2} = 'Start Time'; 
    SvIfo{1, width(SvIfo)+1} = char(t0);
    SvIfo{2, width(SvIfo)-1} = 'Sample Rate (Hz)';
    SvIfo{2, width(SvIfo)} = num2str(SampleRate);
    writecell(SvIfo, svname_+" - info.csv");
    end

    % event table 
    if numel(EV)
    SvEv = table(seconds(EV.Time), seconds(EV.EventEnds), EV.EventLabels, ...
        'VariableNames',{'Start Time (s)', 'End Time (s)', 'Label'});
    writetable(SvEv, svname_+" - events.csv");
    end
end
end

% save figures 
saveasmultiple(fig_stimresponse, svname+" - Stim Response");
saveasmultiple(fig_spec1, svname+" - Spectrum");
for hzn = length(hzns):-1:1
    myfig = fig_auton(hzn);
    if strcmp(class(myfig), 'double')
        myfig = figure(myfig);
    end
    saveasmultiple(myfig, svname+" - Interchannel "+num2str(round(1000*hzns(hzn)))+"ms");
end

%% helpers 

function tblslist = selectTbls(tblsOrig, tblsArtrem, chanselmade, channamesel, chanselind)
% Select between original and art removed table(s). If selected channels
% have not been artifact removed, this will try to pull the original data,
% in which case output table cols might not be in the same order. 
tblsOrig = tblsOrig(:,1);
if ~isempty(tblsArtrem)
    tblslist = tblsArtrem;
    if chanselmade
        for Ti = 1:length(tblslist)
            T = tblslist{Ti};
            for ch = string(T.Properties.VariableNames)
                if ~sum(strcmp(channamesel, ch))
                    % current channel is not desired 
                    T = removevars(T, ch);
                end
            end
            for ch = channamesel
                if ~sum(strcmp(T.Properties.VariableNames, ch))
                    % current channel ch is desired but not art removed
                    terr = nan(size(tblsOrig));
                    for Tj = 1:length(tblsOrig)
                        To = tblsOrig{Tj};
                        terr(Tj) = sum(...
                            abs(seconds(min(T.Time)-min(To.Time)) + ...
                            abs(seconds(max(T.Time)-max(To.Time)))));
                    end
                    [~,Tji] = min(terr); % find the original table with most time overlap
                    To = tblsOrig{Tji};
                    To_ch = array2timetable(To.(ch), ...
                        "VariableNames",ch, "RowTimes",To.Time);
                    T = synchronize(T,To_ch,'first','nearest');
                end
            end
            tblslist{Ti} = T;
        end
    end
else
    tblslist = tblsOrig;
    if chanselmade
        for Ti = 1:length(tblslist)
            T = tblslist{Ti};
            T = T(:,chanselind);
            tblslist{Ti} = T;
        end
    end
end
end

function dcor = DistanceCorrelation(x, y)
% DISTANCE_CORRELATION Compute the distance correlation between x and y.
%
%   dcor = distance_correlation(x, y)
%
%   Inputs:
%       x, y : numeric vectors or matrices with the same number of rows (samples)
%   Output:
%       dcor : scalar, distance correlation coefficient in [0, 1]
%
%   Reference: Szekely, G.J., Rizzo, M.L. (2007).
%   "Measuring and Testing Dependence by Correlation of Distances".
%   Annals of Statistics 35(6):2769–2794.

    % Ensure column vectors
    x = x(:);
    y = y(:);

    n = numel(x);
    if numel(y) ~= n
        error('x and y must have the same length.');
    end

    % limit size for memory and runtime 
    nmax = 2^10;
    if n > nmax
        smp = randperm(n, nmax);
        x = x(smp); y = y(smp);
    end

    % Compute pairwise Euclidean distances
    a = abs(x-x.');
    b = abs(y-y.');

    % Double-center the distance matrices
    A = a - mean(a,1) - mean(a,2) + mean(a(:));
    B = b - mean(b,1) - mean(b,2) + mean(b(:));

    % Compute squared distance covariance and variances
    dcov2 = mean(A(:) .* B(:));
    dvarx2 = mean(A(:) .* A(:));
    dvary2 = mean(B(:) .* B(:));

    % Distance correlation
    dcor = 0;
    if dvarx2 > 0 && dvary2 > 0
        dcor = sqrt(dcov2 / sqrt(dvarx2 * dvary2));
    end
end

function saveasmultiple(fig, filename)
if isvalid(fig)
saveas(fig, filename, 'fig'); % original matlab figure
saveas(fig, filename, 'png'); % preview
else
    warning('figure handle not valid; may have been closed.')
end
end