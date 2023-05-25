function plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
% Creates tf plots for the zmaps with a black outline along the contour of
% the mask_in input

    contourf(x,y,zmap,40,'linecolor','none')
    hold on
    [~,cont] =  contour(x,y,mask_in,1,'linecolor','k')
    cont.LineWidth = 2;
    c = colorbar;
    c.Label.String = 'Z score';
    axis square
    set(gca,'clim',clims_in);
    title(title_in)
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    colormap jet
end


