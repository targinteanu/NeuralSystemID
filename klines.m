function [labels, lines] = klines(X, K, maxIter)

N = size(X,1);
ndim = size(X,2);

% Initialize: random assignments
%labels = randi(K, N, 1);
% Initialize using a GMM 
GM = fitgmdist(X,K,'CovarianceType','full','Replicates',20);
labels = cluster(GM,X);

for iter = 1:maxIter
    
    % Step 1: fit lines
    lines = cell(K,1);
    for k = 1:K
        pts = X(labels==k,:);
        if size(pts,1) < 2
            continue;
        end
        
        mu = mean(pts,1);
        [V,~] = eig(cov(pts));
        dir = V(:,end); % principal direction
        
        lines{k}.point = mu;
        lines{k}.dir = dir / norm(dir);
    end
    
    % Step 2: reassign points
    for i = 1:N
        x = X(i,:)';
        best_k = 1;
        best_dist = inf;
        
        for k = 1:K
            if isempty(lines{k}), continue; end
            
            p = lines{k}.point';
            d = lines{k}.dir;
            
            % orthogonal distance to line
            dist = norm((eye(ndim)-d*d')*(x-p));
            
            if dist < best_dist
                best_dist = dist;
                best_k = k;
            end
        end
        
        labels(i) = best_k;
    end
end