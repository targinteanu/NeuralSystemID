% Match ephys channel names with SEEG contacts identified on imaging. May
% also work with ECoG, etc. 
% 
% Inputs: 
%   chnames: list of channel names from ephys recording 

%% interpret channel names from file
chnames = string(upper(chnames)); % standardize 
chnames = chnames(~contains(chnames, "AINP")); % analog inputs not in head
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

% prep to sort contacts (init) 
labels = nan(height(XYZ),1);
K = length(chnamesu);
lines = nan(K,width(XYZ),2);

%% sort contacts 

sortdone = false; 
while ~sortdone
    figure(figbrain);
    labelsPlots = cell(K,3);
    for label = 1:K
        colr = colorwheel(label/K);
        xyz = XYZ(labels==label,:);
        txt = char(64+label) + string(1:height(xyz));
        labelsPlots{label,1} = text(xyz(:,1),xyz(:,2),xyz(:,3),...
            txt, 'Color',colr, 'FontWeight','bold');
    end
    sortdone = input("Accept? [y/n] ","s");
    sortdone = strcmpi(sortdone, 'y');
    if ~sortdone
    K = input("Number of depth electrodes (default "+string(K)+"): ");
    if isempty(K)
        K = length(chnamesu); % use default if no input
    end
    disp('(Re-)trying klines. Please wait...');
    [labelsNew,linesNew] = klines(XYZ,K,20000,50);
    for l = 1:height(labelsPlots)
        for lp = 1:width(labelsPlots)
            delete(labelsPlots{l,lp});
        end
    end
    figure(figbrain);
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
                "off", 'Color',colr);
        linedir = -linedir;
        labelsPlots{label,3} = quiver3(linecen(1), linecen(2), linecen(3), ...
                linedir(1), linedir(2), linedir(3), ...
                "off", 'Color',colr);
        %labelsPlots{label,4} = text(linecen(1),linecen(2),linecen(3), string(label), 'Color',colr);
    end
    sortdone = input("Accept? [y/n] ","s");
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
    end
end

lblcount = arrayfun(@(l) sum(labels==l), 1:K);
countdiff = length(lblcount) - length(chucount);
if countdiff < 0
    lblcount = [lblcount, zeros(1,-countdiff)];
elseif countdiff > 0
    chucount = [chucount, zeros(1,countdiff)];
end

%% find the best match between lblcount and chucount 

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
labels = arrayfun(@(idx) lblidx(idx), labels);