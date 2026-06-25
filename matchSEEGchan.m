function electbl = matchSEEGchan(chnames, electbl, ...
    XYZacpc, brainmeshLacpc, brainmeshRacpc)
% Match ephys channel names with SEEG contacts identified on imaging. May
% also work with ECoG, etc. 
% 
% Inputs: 
%   chnames: list of channel names from ephys recording 
%   electbl: table of electrode details that has fields x, y, z, AFNI,
%            JuBrain, Brainnetome 
% 
% Optional inputs:
%   XYZacpc: acpc X,Y,Z coordinates (columns=dimensions) in acpc space, or
%            fieldtrip elec_acpc_f struct 
%   brainmesh<L/R>acpc: fieldtrip mesh struct for brain R/L hemisphere in
%            acpc space, or path to .mat file containing "mesh" variable 
% If the above are blank or omitted, mni space will be used instead, with
% coordinates from electbl and fieldtrip template brain. 
% 
% Outputs: 
%   electbl with the field Electrode updated to match chnames 

if nargin < 5
    brainmeshRacpc = [];
    if nargin < 4
        brainmeshLacpc = [];
        if nargin < 3
            XYZacpc = [];
        end
    end
end

%% interpret channel names from file
%chIDs = 1:length(chnames);
chnames = sort(string((chnames))); % standardize 
%chIDs = chIDs(~contains(chnames, "AINP")); 
chnames = chnames(~contains(upper(chnames), "AINP")); % analog inputs not in head
%chIDs = chIDs(~contains(chnames, " BF")); % BFs not visible on CT
chnames = chnames(~contains(upper(chnames), " BF")); % BFs not visible on CT

% remove all digits from each channel name (works on string array)
chnamesu = regexprep(upper(chnames), '-REF', '');
chnamesu = regexprep(chnamesu, '\d+', '');
[chnamesu,ia] = unique(chnamesu);
chucount = diff([sort(ia); length(chnames)+1])';

if ~isempty(chnamesu)

%% view contacts 
if isempty(XYZacpc)
    XYZ = [electbl.x, electbl.y, electbl.z];
else
    XYZ = XYZacpc; % Use provided acpc coordinates if available
    if isstruct(XYZ)
        if isfield(XYZ, 'elecpos')
            XYZ = XYZ.elecpos; % Extract electrode positions if struct
        else
            XYZ = XYZ.chanpos;
        end
    end
    if isempty(brainmeshLacpc) || isempty(brainmeshRacpc)
        warning(['Using template brain in MNI space with provided ' ...
            'electrode coordinates in unknown space.'])
    end
end
coeff = pca(XYZ);
coeffdir = coeff(1:3,1);
coeffdirsgn = sign(mean(XYZ)*coeffdir);
coeffdir = coeffdirsgn*coeffdir;
[az,el] = cart2sph(coeffdir(1),coeffdir(2),coeffdir(3)); % angle for best viewing
az = az*180/pi; el = el*180/pi;

% 3D template brain
global figbrain
figbrain = figure;
if isempty(brainmeshLacpc) || isempty(brainmeshRacpc)
    [ftver, ftpath] = ft_version;
    brainmesh_lh = load([ftpath filesep 'template/anatomy/surface_pial_left.mat'], 'mesh');
    brainmesh_lh = brainmesh_lh.mesh; 
    brainmesh_rh = load([ftpath filesep 'template/anatomy/surface_pial_right.mat'], 'mesh');
    brainmesh_rh = brainmesh_rh.mesh; 
else
    brainmesh_lh = brainmeshLacpc; 
    brainmesh_rh = brainmeshRacpc; 
    if ischar(brainmesh_lh) || isstring(brainmesh_lh)
        brainmesh_lh = ft_read_headshape(brainmesh_lh);
    end
    if ischar(brainmesh_rh) || isstring(brainmesh_rh)
        brainmesh_rh = ft_read_headshape(brainmesh_rh);
    end
    if isempty(XYZacpc)
        warning(['Using provided brain mesh in unknown space with XYZ ' ...
            'coordinates from the table likely in MNI space.'])
    end
end

% show brain mesh 
FaceAlpha = .1; % transparency
FaceColor = .9*[1,1,1];
merh = ft_plot_mesh(brainmesh_rh);
hold on
melh = ft_plot_mesh(brainmesh_lh);
merh.FaceColor = FaceColor;
merh.FaceAlpha = FaceAlpha;
melh.FaceColor = FaceColor;
melh.FaceAlpha = FaceAlpha;
view([az, el]);
material dull;
lighting gouraud;
camlight;

plot3(XYZ(:,1), XYZ(:,2), XYZ(:,3), '.r');

%% sort contacts 

% prep to sort contacts (init) 
labels = nan(height(XYZ),1);
K = length(chnamesu);
lines = nan(K,width(XYZ),2);
sortdone = false; 
labelsPlots = {};

while ~sortdone

    %{
    labelsPlots = plotNewLabels(XYZ, K, labels, lines);
    sortdone = input("Accept? [y/n] ","s");
    sortdone = strcmpi(sortdone, 'y');
    %}

    %if ~sortdone
    K = input("Number of depth electrodes (default "+string(K)+"): ");
    if isempty(K)
        K = length(chnamesu); % use default if no input
    end
    disp('(Re)trying klines. Please wait...');
    [labelsNew,linesNew] = klines(XYZ,K,20000,50);
    for l = 1:height(labelsPlots)
        for lp = 1:width(labelsPlots)
            delete(labelsPlots{l,lp});
        end
    end

    labelsPlots = plotNewLabels(XYZ, K, labelsNew, linesNew);
    sortdone = input("Accept for now? [y/n] ","s");
    sortdone = strcmpi(sortdone, 'y');

    if sortdone
        labels = labelsNew;
    else
        for l = 1:height(labelsPlots)
            for lp = 1:width(labelsPlots)
                delete(labelsPlots{l,lp});
            end
        end
    end
    %end
end

%% manually reassign unnamed labels (1) 

% build a checkerboard UI to reassign points between labels
% Prepare data: for each label, list point indices and coordinates
ptIdxPerLabel = arrayfun(@(l) find(labels==l), 1:K, 'UniformOutput', false);
nPerLabel = cellfun(@numel, ptIdxPerLabel);
maxRows = max(nPerLabel);

% Create figure for checkerboard
figCB = uifigure('Name','Label Reassignment','NumberTitle','off','MenuBar','none', ...
    'Scrollable','on', 'UserData',labels);
t = uitable(figCB, 'Data', cell(maxRows,K), 'ColumnEditable', true(1,K), ...
    'ColumnName', arrayfun(@(l) char(64+l), 1:K, 'UniformOutput', false), ...
    'CellSelectionCallback', @cellSel, 'UserData',[], 'Tag','t');
t.Position(3:4) = [min(1200,150+120*K), min(600,40+30*(maxRows+1))];

% Fill table with point identifiers (show point number and original index)
tbl = strings(maxRows,K);
for c = 1:K
    idxs = ptIdxPerLabel{c};
    for r = 1:numel(idxs)
        tbl(r,c) = "#"+string(r)+" (elec"+string(idxs(r))+")";
    end
end
t.Data = cellstr(tbl);

% Instruction and buttons
uicontrol(figCB,'Style','text','String','Select table cells (one or more) then click "Reassign to Label" to move points to a different label.','HorizontalAlignment','left',...
    'Position',[10, t.Position(2)+t.Position(4)+2, t.Position(3)-20, 20]);
lblPopup = uicontrol(figCB,'Style','popupmenu','String',arrayfun(@(l) char(64+l), 1:K, 'UniformOutput', false),...
    'Position',[10,10,120,22], 'Tag','lblPopup');
btnReassign = uicontrol(figCB,'Style','pushbutton','String','Reassign to Label','Position',[140,10,140,22],...
    'Callback',@(s,e) reassignCallback(s,e, XYZ,K));
btnDone = uicontrol(figCB,'Style','pushbutton','String','Done','Position',[290,10,80,22],...
    'Callback',@(s,e) uiresume(figCB));

uiwait(figCB); % wait until Done pressed

% After UI closed, update labels based on internal state (if any changes left)
if isvalid(figCB)
    labels = figCB.UserData; lines = fitlines(XYZ,K,labels);
    close(figCB);
end

%% find the best match between lblcount and chucount 

lblcount = arrayfun(@(l) sum(labels==l), 1:K);
countdiff = length(lblcount) - length(chucount);
if countdiff < 0
    lblcount = [lblcount, zeros(1,-countdiff)];
    lines = [lines; nan(-countdiff,size(lines,2),size(lines,3))];
elseif countdiff > 0
    chucount = [chucount, zeros(1,countdiff)];
    chnamesu = [chnamesu, strings(1,countdiff)];
end

%{
% Find permutation of lblcount that minimizes sum of absolute differences with chucount
n = length(lblcount);
% For small n use full permutation, otherwise use Hungarian assignment on cost matrix
if n <= 10
    permsIdx = perms(1:n); % each row is a permutation
    bestCost = inf;
    bestPerm = 1:n;
    for p = 1:size(permsIdx,1)
        perm = permsIdx(p,:);
        cost = sum(abs(lblcount(perm) - chucount));
        if cost < bestCost
            bestCost = cost;
            bestPerm = perm;
        end
    end
else
    % Construct cost matrix where cost(i,j) = |lblcount(i) - chucount(j)|
    C = abs(lblcount(:) - chucount(:)');
    % Use munkres (Hungarian). If not available, use assignprob via builtin matchpairs if present.
    if exist('matchpairs','file') == 2
        [pairs, totalCost] = matchpairs(C, Inf);
        % matchpairs returns pairs as [row,col]
        bestPerm = zeros(1,n);
        bestPerm(pairs(:,1)) = pairs(:,2);
        bestCost = totalCost;
    else
        % simple greedy fallback: solve linear assignment via MATLAB's assignment problem using lapjv if available,
        % otherwise use a heuristic: sort both and map sorted indices.
        try
            % attempt to use munkres from File Exchange
            assign = munkres(C);
            bestPerm = assign;
            bestCost = sum(C(sub2ind(size(C), 1:n, bestPerm)));
        catch
            % fallback: sort-based mapping
            [~, sL] = sort(lblcount);
            [~, sC] = sort(chucount);
            bestPerm = zeros(1,n);
            bestPerm(sL) = sC;
            bestCost = sum(abs(lblcount(bestPerm) - chucount));
        end
    end
end

% Apply permutation to lblcount to get best match
lblcount_perm = lblcount(bestPerm);
%}

% may be sufficient to just sort both 
[lblcount, lblidx] = sort(lblcount);
[chucount, chuidx] = sort(chucount);
chnamesu = chnamesu(chuidx);
lines = lines(lblidx,:,:);
labels = arrayfun(@(lbl) find(lblidx==lbl), labels);
%lines = fitlines(XYZ,K,labels); % this should not be needed 
plotNewLabels(XYZ, K, labels, lines);

%% manually reassign named labels (2)

assigndone = false;
while ~assigndone

% create a simple UIFigure with K popupmenus to allow selecting chnamesu for each label
figUI = uifigure('Name','Assign Channel Names','Position',[100 100 300 60+30*K], ...
    'Scrollable','on');
popupHandles = gobjects(K,1); cbxHandles = gobjects(K,1);
lblStrings = cellstr(chnamesu); % options
for k = 1:K
    colr = colorwheel(k/K);
    y = 60 + 30*(K-k);
    uil = uilabel(figUI, 'Position',[10 y 60 22], 'Text', ['Label ',char(64+k)], ...
        'FontColor',colr, 'FontWeight','bold');
    popupHandles(k) = uidropdown(figUI, ...
        'Position',[80 y 200 22], ...
        'Items', lblStrings, ...
        'Value', lblStrings{k});
    
    % add a checkbox to the right of the dropdown to allow selecting this label
    cbxHandles(k) = uicheckbox(figUI, ...
        'Position',[80+200+10 y 60 22], ... % to the right of the dropdown
        'Text','Uncertain', ...
        'Value', false);
end

% add Done button to close and collect selections
btnDone = uibutton(figUI,'push', 'Text','Done', 'Position',[100 10 100 30], ...
    'ButtonPushedFcn', @(btn,event) uiresume(figUI));
uiwait(figUI);

% read selections (if figure still valid)
chnamesu2 = strings(size(chnamesu));
markeduncertain = false(size(chnamesu));
if isvalid(figUI)
    for k = 1:K
        sel = popupHandles(k).Value;
        % assign selected name to all channels in this label
        chnamesu2(k) = string(sel);
        markeduncertain(k) = cbxHandles(k).Value;
    end
    close(figUI);
end

chnamesu = chnamesu2;

%% try assigning named labels 

% get expected dist between channels 
D = [];
for label = 1:K
    xyz = XYZ(labels==label,:);
    linecen = squeeze(lines(label,:,1)); 
    linedir = squeeze(lines(label,:,2));
    d = (xyz-linecen)*linedir';
    d = diff(sort(d));
    D = [D; d];
end
Dmed = median(D); Dstd = std(D);

% initial auto coord-name mapping 
[XYZname, missingname] = matchChElec(K, XYZ,labels,lines, chnames,chnamesu, Dmed,Dstd);
markeduncertaineach = markeduncertain(labels)';

% visualize initial results 
figbrain2 = figure; % fresh figure 

% show brain mesh 
figbrainax(1) = subplot(1,2,1);
FaceAlpha = .1; % transparency
FaceColor = .9*[1,1,1];
merh = ft_plot_mesh(brainmesh_rh);
hold on
melh = ft_plot_mesh(brainmesh_lh);
merh.FaceColor = FaceColor;
merh.FaceAlpha = FaceAlpha;
melh.FaceColor = FaceColor;
melh.FaceAlpha = FaceAlpha;
view([az, el]);
material dull;
lighting gouraud;
camlight;

% show labeled elecs 
for label = 1:K
    colr = colorwheel(label/K);
    xyz = XYZ(labels==label,:);
    chnamenumk = XYZname(labels==label);
    plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.', 'Color',colr);
    text(xyz(:,1),xyz(:,2),xyz(:,3), chnamenumk, 'Color',colr, 'FontWeight','bold');
end

missingListStr = strjoin(missingname, ', ');
missingListStr = char(missingListStr); 
missingListStr = missingListStr(3:end); % remove leading comma 
missingListStr = ['Missing: ',missingListStr];
disp(missingListStr);

% show brain mesh 
figbrainax(2) = subplot(1,2,2);
FaceAlpha = .1; % transparency
FaceColor = .9*[1,1,1];
merh = ft_plot_mesh(brainmesh_rh);
hold on
melh = ft_plot_mesh(brainmesh_lh);
merh.FaceColor = FaceColor;
merh.FaceAlpha = FaceAlpha;
melh.FaceColor = FaceColor;
melh.FaceAlpha = FaceAlpha;
view([az, el]);
material dull;
lighting gouraud;
camlight;

% show labeled elecs 
for label = 1:K
    colr = colorwheel(label/K);
    xyz = XYZ(labels==label,:);
    linecen = squeeze(lines(label,:,1)); 
    chnamek = chnamesu(label);
    plot3(xyz(:,1),xyz(:,2),xyz(:,3), ...
        '.', 'Color',colr);
    text(linecen(1),linecen(2),linecen(3), ...
        chnamek, 'Color',colr, 'FontWeight','bold');
    %{
    linelen = [max(xyz(:,1)), max(xyz(:,2)), max(xyz(:,3))] - ...
              [min(xyz(:,1)), min(xyz(:,2)), min(xyz(:,3))];
    linelen = norm(linelen);
    linedir = squeeze(lines(label,:,2)); linedir = linedir*linelen/2;
    quiver3(linecen(1), linecen(2), linecen(3), ...
            linedir(1), linedir(2), linedir(3), ...
            "off", 'Color',colr, 'LineWidth',1);
    linedir = -linedir;
    quiver3(linecen(1), linecen(2), linecen(3), ...
            linedir(1), linedir(2), linedir(3), ...
            "off", 'Color',colr, 'LineWidth',1);
    %}
end

linkprop(figbrainax, {'CameraPosition', 'CameraTarget', 'CameraUpVector'});

disp([XYZname, electbl.AFNI, electbl.JuBrain, electbl.Brainnetome])

    assigndone = input("Accept final? [y/n] ","s");
    assigndone = strcmpi(assigndone, 'y');

end

electbl.Electrode = XYZname; 
electbl.MatchUncertain = markeduncertaineach;

end
end

%% helper(s) 

function cnum = getnumfromstr(cname)
cname = char(cname); 
cnum = cname((cname >= 48) & (cname <= 57)); 
%cnum = str2double(cnum);
end

% re-fit lines from existing points without re-running klines 
function lines = fitlines(X,K,labels)
N = size(X,1);
ndim = size(X,2);
    % Step 1: fit lines
    lines = nan(K,ndim,2);
    for k = 1:K
        pts = X(labels==k,:);
        if size(pts,1) < 2
            continue;
        end
        
        mu = mean(pts,1);
        [V,~] = eig(cov(pts));
        dir = V(:,end); % principal direction
        
        lines(k,:,1) = mu;
        lines(k,:,2) = dir' / norm(dir);
    end
end

function [XYZname, missingname] = matchChElec(K, XYZ,labels,lines, chnames,chnamesu, Dmed,Dstd)
%XYZrow = (1:size(XYZ,1))'; 
XYZname = strings(size(XYZ,1),1);
missingname = "";
for label = 1:K

    % get names/IDs from ephys 
    chnamek = chnamesu(label);
    chsel = contains(chnames,chnamek);
    chnamenumk = chnames(chsel); %chIDk = chIDs(chsel); 
    chnumk = arrayfun(@(str) getnumfromstr(str), chnamenumk, 'UniformOutput',false);
    chnumk = string(chnumk);
    % remove names that contain other names (e.g. LA vs LAH)
    chnamesk = chnamenumk; 
    for ch = 1:length(chnamesk)
        chnamesk_ch = split(chnamenumk(ch), chnumk(ch));
        chnamesk(ch) = chnamesk_ch(1);
    end
    chnumk = str2double(chnumk);
    chsel = strcmp(chnamesk, chnamek); 
    chnamenumk = chnamenumk(chsel); %chIDk = chIDk(chsel); 
    chnumk = chnumk(chsel); chnamesk = chnamesk(chsel);
    % sort and handle duplicates
    idxtf = false(size(chnumk));
    [chnumk,idx] = unique(chnumk); % treat duplicates as missing
    idxtf(idx) = true;
    missingname = [missingname; chnamenumk(~idxtf)'];
    chnamenumk = chnamenumk(idx); %chIDk = chIDk(idx);
    [chnumk,idx] = sort(chnumk, 'ascend'); % assume numbering is spatial; lowest=deepest
    chnamenumk = chnamenumk(idx); %chIDk = chIDk(idx);

    % get ordered positions 
    xyz = XYZ(labels==label,:); %r = XYZrow(labels==label,:);
    linecen = squeeze(lines(label,:,1)); 
    linedir = squeeze(lines(label,:,2));
    rightsided = linecen(1) > 0;
    rightpointing = linedir(1) > 0;
    outpointing = rightsided == rightpointing;
    if ~outpointing
        linedir = -linedir; % force outpointing
    end
    d = (xyz-linecen)*linedir'; % lower (more negative) = deeper
    [d,didx] = sort(d); % order wrt same label elecs
    r = (xyz-mean(XYZ))./std(XYZ); % diff from all elecs 
    r = rms(r,2);

    % any channels not seen on imaging?  
    dd = diff(d); ddout = (dd-Dmed) > 3*Dstd;
    while length(chnamenumk) > length(d)
        if any(ddout)
            [~,ddi] = max(dd);
            % assume the missing chan was between ddi+1 and ddi
            missingname = [missingname; chnamenumk(ddi+1)];
            chnamenumk = chnamenumk([1:ddi, (ddi+2):end]);
            chnumk = chnumk([1:ddi, (ddi+2):end]);
            %chIDk = chIDk([1:ddi, (ddi+2):end]);
            % plug ddi so the same index doesn't get flagged again 
            dd(ddi) = -inf; ddout(ddi) = false;
        else
            % assume the missing chan is most superficial
            missingname = [missingname; chnamenumk(end)];
            chnamenumk = chnamenumk(1:(end-1));
            chnumk = chnumk(1:(end-1));
            %chIDk = chIDk(1:(end-1));
        end
    end

    % any imaged "channels" not recorded? 
    dd = diff(d); ddout = (Dmed-dd) > 3*Dstd;
    while length(d) > length(chnamenumk)
        if any(ddout)
            [~,ddi] = min(dd);
            % ddi+1 is too close to ddi to be distinct 
            d = d([1:ddi, (ddi+2):end]); 
            didx = didx([1:ddi, (ddi+2):end]);
            dd = diff(d); ddout = (Dmed-dd) > 3*Dstd;
        else
            % assume superficial-most contact is incorrect 
            chnamenumk = [chnamenumk,"?"];
            chnumk = [chnumk,nan]; %chIDk = [chIDk,nan];
            %{
            [~,rj] = min(r); % try to put this near the middle 
            % mark most eccentric elec as unknown 
            [~,ri] = max(r); rij = ri-rj;
            sj = ceil(length(chnamenumk)/2); si = sj+rij;
            if si < 1
                chnamenumk = ["?",chnamenumk];
                chnumk = [nan,chnumk]; %chIDk = [nan,chIDk];
            elseif si > length(chnamenumk)
                chnamenumk = [chnamenumk,"?"];
                chnumk = [chnumk,nan]; %chIDk = [chIDk,nan];
            else
                chnamenumk = [chnamenumk(1:(si-1)),"?",chnamenumk((si):end)];
                chnumk = [chnumk(1:(si-1)),nan,chnumk((si):end)];
                %chIDk = [chIDk(1:(si-1)),nan,chIDk((si):end)];
            end
            % plug ri so the same index doesn't get flagged again 
            r(ri) = nan;
            %}
        end
    end

    % assign pairs 
    % everything is now the side of d 
    labelsk = find(labels==label);
    labelsk = labelsk(didx);
    XYZname(labelsk) = chnamenumk';
end
end


% update the labels/lines/text on figbrain
function labelsPlots = plotNewLabels(XYZ, K, labelsNew, linesNew)
    global figbrain
    figure(figbrain);
    figgraphics = gca().Children;
    for gr = figgraphics'
        if isprop(gr, 'Type')
            grType = gr.Type;
            if strcmpi(grType,'Text') || strcmpi(grType,'Quiver')
                delete(gr);
            end
        end
    end
    labelsPlots = cell(K,3);
    for label = 1:K
        colr = colorwheel(label/K);
        xyz = XYZ(labelsNew==label,:);
        txt = char(64+label) + string(1:height(xyz));
        labelsPlots{label,1} = text(xyz(:,1),xyz(:,2),xyz(:,3),...
            txt, 'Color',colr, 'FontWeight','bold');
        linelen = [max(xyz(:,1)), max(xyz(:,2)), max(xyz(:,3))] - ...
                  [min(xyz(:,1)), min(xyz(:,2)), min(xyz(:,3))];
        linelen = norm(linelen);
        linecen = squeeze(linesNew(label,:,1));
        linedir = squeeze(linesNew(label,:,2)); linedir = linedir*linelen/2;
        labelsPlots{label,2} = quiver3(linecen(1), linecen(2), linecen(3), ...
                linedir(1), linedir(2), linedir(3), ...
                "off", 'Color',colr, 'LineWidth',1);
        linedir = -linedir;
        labelsPlots{label,3} = quiver3(linecen(1), linecen(2), linecen(3), ...
                linedir(1), linedir(2), linedir(3), ...
                "off", 'Color',colr, 'LineWidth',1);
        %labelsPlots{label,4} = text(linecen(1),linecen(2),linecen(3), string(label), 'Color',colr);
    end
end


% Callback: record selected cells
function cellSel(src, event, selCells)
    if isempty(event.Indices)
        selCells = [];
    else
        selCells = unique(event.Indices, 'rows');
    end
    src.UserData = selCells;
end

% Callback: perform reassignment of selected points to chosen label
function reassignCallback(src, evt, XYZ,K)
    global figbrain
    labels = src.Parent.UserData;
    %lblPopup = src.Parent.Children(3); % change to tag?
    %t = src.Parent.Children(end); % change to tag?
    figgraphics = src.Parent.Children;
    figgraphicsnames = {figgraphics.Tag};
    lblPopup = figgraphics(strcmp(figgraphicsnames, 'lblPopup'));
    t = figgraphics(strcmp(figgraphicsnames, 't'));
    selCells = t.UserData;
    if isempty(selCells)
        uialert(src.Parent,'No cells selected.','Warning','Icon','warning');
        return;
    end
    newLabel = lblPopup.Value; % numeric label target
    % Determine point indices from selected cells: parse strings like "pt3 (#12)"
    data = srcTableData(t);
    moved = false;
    for k = 1:size(selCells,1)
        r = selCells(k,1); c = selCells(k,2);
        cellStr = string(data{r,c});
        if strlength(cellStr)==0
            continue;
        end
        % extract number inside parentheses
        tok = sscanf(cellStr, '#%u (elec%u)');
        if isempty(tok), continue; end
        ptIdx = tok(2);
        if isnan(ptIdx), continue; end
        % reassign label array
        labels(ptIdx) = newLabel;
        moved = true;
    end
    if moved
        % rebuild table contents to reflect new grouping
        ptIdxPerLabel = arrayfun(@(l) find(labels==l), 1:K, 'UniformOutput', false);
        nPerLabel = cellfun(@numel, ptIdxPerLabel);
        maxRows = max(nPerLabel);
        tbl = strings(maxRows,K);
        for cc = 1:K
            idxs = ptIdxPerLabel{cc};
            for rr = 1:numel(idxs)
                tbl(rr,cc) = "#"+string(rr)+" (elec"+string(idxs(rr))+")";
            end
        end
        t.Data = cellstr(tbl);
        selCells = []; % clear selection
        t.UserData = selCells;
        % update labels list 
        src.Parent.UserData = labels;
        % rebuild lines 
        lines = fitlines(XYZ,K,labels);
        % Also update plot labels on brain figure
        figure(figbrain);
        plotNewLabels(XYZ,K,labels,lines); 
    end
end

% helper to safely get table data as cell array of strings
function c = srcTableData(tblHandle)
    d = tblHandle.Data;
    if iscell(d)
        c = d;
    else
        % table.Data might be string array; convert to cellstr with possibly empty cells
        c = cell(size(d));
        for ii = 1:numel(d)
            c{ii} = char(d(ii));
        end
        c = reshape(c, size(d));
    end
end