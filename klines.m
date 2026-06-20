function [labels, lines] = klines(X, K, maxIter, gmReplicates)

N = size(X,1);
ndim = size(X,2);

% Initialize: random assignments
%labels = randi(K, N, 1);
% Initialize using a GMM 
original_state = warning('off','all');
GM = fitgmdist(X,K,'CovarianceType','full','Replicates',gmReplicates);
labels = cluster(GM,X);
warning(original_state);

for iter = 1:maxIter
    
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
    
    % Step 2: reassign points
    for i = 1:N
        x = X(i,:)';
        best_k = 1;
        best_dist = inf;
        
        for k = 1:K
            if any(isnan(lines(k,:,:)),'all'), continue; end
            
            p = lines(k,:,1)';
            d = lines(k,:,2)';
            
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