function learnrate = LearnrateEst(dta, Lmax)
% 
% Determine a learning rate (step size) for gradient descent on input data
% dta. Learning rate is estimated using the autocorrelation of subsamples of
% dta no greater than Lmax. dta must be a column of data. With Lmax, longer
% values will get better results but take a long time
% 

%% handle inputs 
if nargin < 2
    Lmax = [];
end
if isempty(Lmax)
    Lmax = 5000;
end
if Lmax < 3
    error('Lmax too small.')
end

%% first subsampling & eigenvalue estimate 
x = dta; L = Lmax;
if length(x) > L
    % subsample data
    x = sort(x); % get as many different vals as possible. eig should be invariant
    xi = linspace(1,length(x),L);
    xi = round(xi); xi = unique(xi);
    x = x(xi);
end
l1 = eig(x*x'); L1 = L;

%% second subsampling & eigenvalue estimate 
x = dta; L = ceil(.5*Lmax);
if length(x) > L
    % subsample data
    x = sort(x); % get as many different vals as possible. eig should be invariant
    xi = linspace(1,length(x),L);
    xi = round(xi); xi = unique(xi);
    x = x(xi);
end
l2 = eig(x*x'); L2 = L;

%% max eigenvalue & learning rate calculation

% assume max eigenvalue is a linear function of sample size 
l1max = max(l1); l2max = max(l2); 
m = (l1max-l2max)/(L1-L2);
lmax = l1max + m*(length(dta)-L1);

% choose learn rate based on max eigenvalue 
learnrate = 2/lmax; % max stable value 

end