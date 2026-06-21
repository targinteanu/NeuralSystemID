% Match ephys channel names with SEEG contacts identified on imaging. May
% also work with ECoG, etc. 
% 
% Inputs: 
%   chnames: list of channel names from ephys recording 

%% interpret channel names from file
chIDs = 1:length(chnames);
chnames = string(upper(chnames)); % standardize 
chIDs = chIDs(~contains(chnames, "AINP")); 
chnames = chnames(~contains(chnames, "AINP")); % analog inputs not in head
chIDs = chIDs(~contains(chnames, " BF")); % BFs not visible on CT
chnames = chnames(~contains(chnames, " BF")); % BFs not visible on CT

% remove all digits from each channel name (works on string array)
chnamesu = regexprep(chnames, '-REF', '');
chnamesu = regexprep(chnamesu, '\d+', '');
[chnamesu,ia] = unique(chnamesu);
chucount = diff([ia; length(chnames)+1])';

%% view contacts 
XYZ = [electbl.x, electbl.y, electbl.z];
coeff = pca(XYZ);
[az,el] = cart2sph(coeff(1,1),coeff(2,1),coeff(3,1)); % angle for best viewing
az = az*180/pi; el = el*180/pi;

% 3D template brain
global figbrain
figbrain = figure;
[ftver, ftpath] = ft_version;
load([ftpath filesep 'template/anatomy/surface_pial_left.mat']);
template_lh = mesh; clear mesh;
load([ftpath filesep 'template/anatomy/surface_pial_right.mat']);
template_rh = mesh; clear mesh;

% show brain mesh 
FaceAlpha = .1; % transparency
FaceColor = .9*[1,1,1];
merh = ft_plot_mesh(template_rh);
hold on
melh = ft_plot_mesh(template_lh);
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
elseif countdiff > 0
    chucount = [chucount, zeros(1,countdiff)];
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
%plotNewLabels(XYZ, K, labels, lines);

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

% initial coord-name mapping 
%XYZrow = (1:size(XYZ,1))'; 
XYZname = strings(size(XYZ,1),1);
missingname = "";
for label = 1:K

    % get names/IDs from ephys 
    chnamek = chnamesu(label);
    chIDk = chIDs(contains(chnames,chnamek));
    chnamenumk = chnames(contains(chnames,chnamek));
    chnumk = arrayfun(@(str) getnumfromstr(str), chnamenumk);
    idxtf = false(size(chnumk));
    [chnumk,idx] = unique(chnumk); % treat duplicates as missing
    idxtf(idx) = true;
    missingname = [missingname; chnamenumk(~idxtf)'];
    chnamenumk = chnamenumk(idx); chIDk = chIDk(idx);
    [chnumk,idx] = sort(chnumk, 'ascend'); % assume numbering is spatial; lowest=deepest
    chnamenumk = chnamenumk(idx); chIDk = chIDk(idx);

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
            chIDk = chIDk([1:ddi, (ddi+2):end]);
            % plug ddi so the same index doesn't get flagged again 
            dd(ddi) = -inf; ddout(ddi) = false;
        else
            % assume the missing chan is most superficial
            missingname = [missingname; chnamenumk(end)];
            chnamenumk = chnamenumk(1:(end-1));
            chnumk = chnumk(1:(end-1));
            chIDk = chIDk(1:(end-1));
        end
    end

    % any imaged channels not recorded? 
    dd = diff(d); ddout = (Dmed-dd) > 3*Dstd;
    while length(d) > length(chnamenumk)
        if any(ddout)
            [~,ddi] = min(dd);
            % ddi+1 is too close to ddi to be distinct 
            d = d([1:ddi, (ddi+2):end]); 
            didx = didx([1:ddi, (ddi+2):end]);
            dd = diff(d); ddout = (Dmed-dd) > 3*Dstd;
        else
            [~,rj] = min(r); % try to put this near the middle 
            % mark most eccentric elec as unknown 
            [~,ri] = max(r); rij = ri-rj;
            sj = ceil(length(chnamenumk)/2); si = sj+rij;
            if si < 1
                chnamenumk = ["?",chnamenumk];
                chnumk = [nan,chnumk]; chIDk = [nan,chIDk];
            elseif si > length(chnamenumk)
                chnamenumk = [chnamenumk,"?"];
                chnumk = [chnumk,nan]; chIDk = [chIDk,nan];
            else
                chnamenumk = [chnamenumk(1:(si-1)),"?",chnamenumk((si):end)];
                chnumk = [chnumk(1:(si-1)),nan,chnumk((si):end)];
                chIDk = [chIDk(1:(si-1)),nan,chIDk((si):end)];
            end
            % plug ri so the same index doesn't get flagged again 
            r(ri) = nan;
        end
    end

    % assign pairs 
    % everything is now the side of d 
    labelsk = find(labels==label);
    labelsk = labelsk(didx);
    XYZname(labelsk) = chnamenumk';
end

% visualize initial results 
figbrain2 = figure; % fresh figure 

% show brain mesh 
FaceAlpha = .1; % transparency
FaceColor = .9*[1,1,1];
merh = ft_plot_mesh(template_rh);
hold on
melh = ft_plot_mesh(template_lh);
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

%% manually reassign named labels (2)


%% helper(s) 

function cnum = getnumfromstr(cname)
cname = char(cname); 
cnum = cname((cname >= 48) & (cname <= 57)); 
cnum = str2double(cnum);
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