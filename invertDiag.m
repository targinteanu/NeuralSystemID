function [Minv, Minv_cell] = invertDiag(M_cell)
% Input a block-diagonal matrix M as a cell array where each cell is a
% block; return M inverted, both in the same cell form and in matrix form. 

Minv_cell = M_cell; 
for r = 1:length(Minv_cell)
    Mr = Minv_cell{r};
    if (height(Mr)==1) && (width(Mr)==1)
        Mr_inv = 1/Mr;
    elseif (height(Mr)==2) && (width(Mr)==2)
        a = Mr(1,1); b = Mr(1,2); c = Mr(2,1); d = Mr(2,2);
        Mr_inv = (1/(a*d - b*c)) * [d, -b; -c, a];
    else
        Mr_inv = Mr^-1;
    end
    Minv_cell{r} = Mr_inv;
end

h = cellfun(@height, Minv_cell); H = sum(h);
Minv = zeros(H); 
Minv = Minv * M_cell{1}(1); % converts to correct type
r = 1; 
for c = 1:length(Minv_cell)
    Mr_inv = Minv_cell{c};
    R = height(Mr_inv)-1;
    Minv(r:(r+R), r:(r+R)) = Mr_inv;
    r = r+R+1;
end

end