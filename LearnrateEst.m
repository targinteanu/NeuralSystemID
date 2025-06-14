function learnrate = LearnrateEst(x, N, Lmax)
%
% Determine a learning rate (step size) for gradient descent on input data
% x with N taps. Learning rate is estimated using the autocorrelation of
% subsamples of x no greater than Lmax. x must be a column of data. With
% Lmax, longer values will get better results but take a long time
%

% handle inputs 
if nargin < 3
    Lmax = [];
end
if isempty(Lmax)
    Lmax = 1e6;
end
Lmax = min(Lmax, length(x));

% estimate autocorrelation
RR = nan(N,N,Lmax);
for t = N:Lmax
    xt = x((t-N+1):t);
    R = xt*xt';
    RR(:,:,t) = R;
end
figure;
R = mean(RR,3,'omitnan'); 
l = eig(R); lmax = max(l);

% determine learn rate uring autocorr
learnrate = 2/lmax; % max stable value 
learnrate = .5*learnrate; % for safety

end