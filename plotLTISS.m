function [fig, tl] = plotLTISS(sys)

A = sys.A; B = sys.B; C = sys.C; D = sys.D;
W1 = width(A); W2 = width(B); H1 = height(A); H2 = height(C);

% force min ratios to 1:5
N = 5; 
if W2/(W1+W2) < 1/N
    W2 = 1; W1 = N-1;
elseif W1/(W1+W2) < 1/N
    W1 = 1; W2 = N-1;
end
if H2/(H1+H2) < 1/N
    H2 = 1; H1 = N-1;
elseif H1/(H1+H2) < 1/N
    H1 = 1; H2 = N-1;
end

fig = figure; 
tl = tiledlayout(H1+H2, W1+W2);

nexttile([H1, W1]); imagesc(A); title('A'); colorbar;
nexttile([H1, W2]); imagesc(B); title('B'); colorbar;
nexttile([H2, W1]); imagesc(C); title('C'); colorbar;
nexttile([H2, W2]); imagesc(D); title('D'); colorbar;

end