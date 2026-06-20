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
figure;
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

labels = nan(height(XYZ),1);
sortdone = false; 
K = length(chnamesu);
while ~sortdone
    labelsPlots = cell(K,1);
    for label = 1:K
        xyz = XYZ(labels==label,:);
        labelsPlots{label} = plot3(xyz(:,1),xyz(:,2),xyz(:,3),...
            'o','Color',colorwheel(label/K));
    end
    K = input("Number of depth electrodes (default "+string(K)+"): ");
    if isempty(K)
        K = length(chnamesu); % use default if no input
    end
    labelsNew = klines(XYZ,K,10000);
    for l = 1:length(labelsPlots)
        delete(labelsPlots{l});
    end
    labelsPlots = cell(K,1);
    for label = 1:K
        xyz = XYZ(labelsNew==label,:);
        labelsPlots{label} = plot3(xyz(:,1),xyz(:,2),xyz(:,3),...
            'o','Color',colorwheel(label/K));
    end
    sortdone = input("Accept? [y/n] ","s");
    sortdone = strcmpi(sortdone, 'y');
    if sortdone
        labels = labelsNew;
    else
        for l = 1:length(labelsPlots)
            delete(labelsPlots{l});
        end
    end
end

lblcount = arrayfun(@(l) sum(labels==l), 1:K);