function [fig, tl] = plotLTISS(sys)

A = sys.A; B = sys.B; C = sys.C; D = sys.D;

fig = figure; 
tl = tiledlayout(height(A)+height(C), width(A)+width(B));

nexttile([height(A), width(A)]); imagesc(A); colorbar;
nexttile([height(B), width(B)]); imagesc(B); colorbar;
nexttile([height(C), width(C)]); imagesc(C); colorbar;
nexttile([height(D), width(D)]); imagesc(D); colorbar;

end