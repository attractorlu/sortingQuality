function plotAsProbe(chanAmps, xc, yc)

cmap = colormap('cool');
ncmap = size(cmap,1);

amp = chanAmps;
amp = amp-min(amp);
amp = amp/max(amp);
color_ind = round(amp * (ncmap-1) ) + 1;

scatter(xc, yc, 20, cmap(color_ind,:), 'filled')

% xlim( [min(xc) max(xc)] * 1.2 )
% ylim( [-100 max(yc) * 1.1] )

