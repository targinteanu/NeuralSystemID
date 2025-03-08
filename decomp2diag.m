function [V, D, Vcell, Dcell] = decomp2diag(M)
% Decompose a matrix M to V * D * V^-1 where D is block diagonal and real. 
% If M has a full set of eigenvalues where complex eigenvalues are paired
% with complementary eigenvectors, D will have blocks of size 2 or lower. 
% Dcell is a cell array of ordered blocks. 

tol = 10*eps; % tolerance for magnitude to consider = 0

[W,E] = eig(M); E = diag(E)';
Er = real(E); Ei = imag(E);
EisReal = abs(imag(E)) <= tol;

% handle real eigenvalues (no processing) 
VReal = W(:,EisReal); DReal = E(EisReal);
VcellReal = mat2cell(VReal,height(VReal),ones(1,width(VReal)));
DcellReal = num2cell(DReal);
% check for remaining complex values 
if sum(abs(imag(VReal)) > tol)
    warning('M has complex eigenvectors with real eigenvalues!')
else
    VReal = real(VReal);
end

% handle complex eigenvalues
W = W(:,~EisReal); E = E(~EisReal);
VComp = []; 
VcellComp = {}; DcellComp = {};
while ~isempty(E)
    e = E(1); % current eigenvalue 
    w = W(:,1); % current eigenvector 
    E = E(2:end); W = W(:,2:end);
    e_ = conj(e); 
    [~,e_i] = min(abs(E - e_));
    e_ = E(e_i); w_ = W(:,e_i); % complex conjugates 
    if abs(e - conj(e_)) > tol
        error('Some complex eigenvalues do not have conjugates.')
    end
    E = E([1:(e_i-1), (e_i+1):end]); W = W(:,[1:(e_i-1), (e_i+1):end]);
    % if vals are a +- bi and vecs are u +- iv: 
    % D should get a block of [a, b; -b, a] 
    % and V should get a block of [u, v]
    u = w + w_; u = .5*real(u);
    v = w - w_; v = .5*imag(v);
    a = e + e_; a = .5*real(a);
    b = e - e_; b = .5*imag(b);
    VV = [u, v]; DD = [a, -b; b, a];
    VComp = [VComp, VV]; VcellComp = [VcellComp, VV]; 
    DcellComp = [DcellComp, DD]; 
end

% output 
V = [VReal, VComp];  
Vcell = [VcellReal, VcellComp]; Dcell = [DcellReal, DcellComp];
D = zeros(size(M));
r = 1;
for c = 1:length(Dcell)
    DD = Dcell{c};
    R = height(DD)-1;
    D(r:(r+R), r:(r+R)) = DD;
    r = r+R+1;
end

end