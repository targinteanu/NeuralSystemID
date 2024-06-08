function [srcness, snkness, fig] = SourceSink(A, chanlocs)

srcness = sqrt(sum(A.^2,1)); 
snkness = sqrt(sum(A.^2,2));

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

fig(2) = figure;
for t = 1:T
    subplot(2,T,t);
    topoplot(srcness(:,:,t), chanlocs, 'maplimits', 'maxmin'); colorbar;
    subplot(2,T,t+T);
    topoplot(snkness(:,:,t), chanlocs, 'maplimits', 'maxmin'); colorbar; 
end

end