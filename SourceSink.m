function [srcness, snkness, fig] = SourceSink(A, chanlocs)

srcness = sqrt(sum(A.^2,1)); 
snkness = sqrt(sum(A.^2,2));
%srcrng = [min(srcness(:)), max(srcness(:))]; 
%snkrng = [min(snkness(:)), max(snkness(:))]; 
srcrng = 'maxmin'; snkrng = 'maxmin';

chlbl = {chanlocs.labels};
T = size(A,3);
fig(1) = figure;
for t = 1:T
    text(srcness(:,:,t)/sum(srcness(:,:,t)), ...
         snkness(:,:,t)/sum(snkness(:,:,t)), ...
         chlbl, ...
         'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
         'Color',colorwheel(t/T));
    hold on; 
end
grid on; 
ylabel('Sink-ness'); xlabel('Source-ness');

fig(2) = figure('Units','normalized', 'Position',[.05,.05,.9,.9]);
if T <= 8
    for t = 1:T
        subplot(2,T,t);
        topoplot(srcness(:,:,t), chanlocs, 'maplimits', srcrng); colorbar;
        subplot(2,T,t+T);
        topoplot(snkness(:,:,t), chanlocs, 'maplimits', snkrng); colorbar; 
    end
else
    W = ceil(sqrt(2*T)); 
    if mod(W,2)
        W = W+1;
    end
    H = ceil(T/W); TT = H*W;
    for t = 1:T
        subplot(2*H,W,t); 
        topoplot(srcness(:,:,t), chanlocs, 'maplimits', srcrng); colorbar;
        subplot(2*H,W,t+TT);
        topoplot(snkness(:,:,t), chanlocs, 'maplimits', snkrng); colorbar;
    end
end

end