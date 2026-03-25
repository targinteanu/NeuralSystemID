function [score,coeff,coeffInv,Xr,mu] = PCAmanifold(X, cutoff)
% 
% Find a lower-dim manifold using PCA. 
% Input data X as column matrix or (time)table. Optionally, input a cutoff
% value for % variance explained by manifold. 
% Outputs: 
%   score: PC scores of selected components on manifold
%   coeff: selected component weightings such that (X-mu)*coeff = score
%   coeffInv: can be used to reconstruct X ~ (score*coeffInv)+mu
%   Xr: original data reconstructed from selected components. 
%   mu: original data mean values 

if nargin < 2
    cutoff = [];
end

t = [];
if istimetable(X)
    t = X.Time; 
    X = table2array(X);
elseif istable(X)
    X = table2array(X);
end

H = height(X); W = width(X); 
if W > H
    error('Cannot have more variables than samples.')
end

% detect outliers 
OLmask = false(H,1); % replace with mask of outlier samples to ignore 

% run PCA 
[coeff,score,latent,tsquared,explained,mu] = pca(X(~OLmask,:));

% if cutoff specified, treat it as a threshold for cumulative explained
Expl = cumsum(explained);
if ~isempty(cutoff)
    q = find(Expl >= cutoff, 1);
    if isempty(q)
        q = W;
    end
else
    q = sum(latent > (H-1)*eps(latent(1)));
    % this is the method used by pca. Alternatively, use rank(X,tol)
end
disp("Processing "+string(q)+" of "+string(W)+" components, explaining "...
    +string(Expl(q))+"%.");

score = score(:,1:q); 
coeffInv = coeff^-1; coeffInv = coeffInv(1:q,:);
coeff = coeff(:,1:q);
Xr = (score*coeffInv)+mu;

end