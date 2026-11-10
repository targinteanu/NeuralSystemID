function [A, k1, c2] = pinkcorrect(A,f)
% correct for noise that obeys Anoise = k1*f^c2
% i.e. ln(Anoise) = c2*ln(f) + c2*ln(k1)
if f(1) < 2*eps
    f0 = 0; f = f(2:end);
    A0 = A(1,:); A = A(2:end,:);
else
    f0 = zeros(0,width(f)); A0 = zeros(0,width(A));
end
lnA = log(A); lnf = log(f); F = [ones(size(lnf)), lnf];
c = F\lnA; 
% c1 = c2*ln(k1), i.e. k1 = exp(c1/c2)
c2 = c(2); k1 = exp(c(1)/c(2));
lnAnoise = F*c;
lnA = lnA - lnAnoise; A = exp(lnA);
A = [A0; A];
end