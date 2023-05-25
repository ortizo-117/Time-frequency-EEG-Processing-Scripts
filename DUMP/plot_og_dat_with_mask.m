function plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in)
% Creates tf plots for the zmaps with a black outline along the contour of
% the mask_in input

    contourf(x,y,data_in,40,'linecolor','none')
    hold on
    [~,cont] =  contour(x1,y,logical(mask_in),1,'linecolor','k')
    cont.LineWidth = 2;
    c = colorbar;
    c.Label.String = '% change from baseline';
    axis square
    set(gca,'clim',clims_in);
    title(title_in)
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    xlim(my_xlims);
    colormap jet
end

