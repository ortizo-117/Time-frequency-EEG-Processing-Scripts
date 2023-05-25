%% cheps vs leps figure maker


%clear all; close all
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\scripts_functions\Colormaps
load('new_color_map_2.mat')
load('custom4.mat')
gamma_in = load('nppc_gamma.mat');
lf_in = load('nppc_lf.mat');
load('gamma_results.mat');
load('lf_results.mat');
all_data = concatenate_all_frequencies(gamma_in,lf_in);
my_results = concatenate_all_results(gamma_results,lf_results);
%%
my_custom_map = plasma(250); %inferno fake parula plasma viridis 
%my_custom_map = customcolormap_preset('brown-white-pool');
x = all_data.x;
y = all_data.y;

clims_in = [-5 5];
f = figure;
hold on
ax(1) = subplot(3,3,1)
title_in = ["C3-Mast","CHEPs"];
zmap = all_data.c1_zmap_c3;
mask_in  = all_data.c1_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,logical(mask_in),x,y,clims_in,title_in)
hold on
ax(2) =subplot(3,3,2)
title_in = ["Cz-Mast","CHEPs"];
zmap = all_data.c1_zmap_cz;
mask_in  = all_data.c1_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,logical(mask_in),x,y,clims_in,title_in)
ax(3) =subplot(3,3,3)
title_in = ["C4-Mast","CHEPs"];
zmap = all_data.c1_zmap_c4;
mask_in  = all_data.c1_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,logical(mask_in),x,y,clims_in,title_in)

ax(4) =subplot(3,3,4)
title_in = ["LEPs"];
zmap = all_data.c3_zmap_c3;
mask_in  = all_data.c3_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,logical(mask_in),x,y,clims_in,title_in)
hold on
ax(5) =subplot(3,3,5)
zmap = all_data.c3_zmap_cz;
mask_in  = all_data.c3_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,logical(mask_in),x,y,clims_in,title_in)
ax(6) =subplot(3,3,6)
zmap = all_data.c3_zmap_c4;
mask_in  = all_data.c3_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,logical(mask_in),x,y,clims_in,title_in)
sgtitle('Frequency response with cluster correction')


my_clims = [-5 5];
my_xlims = [-150 800];
%C3
ax(7) =subplot(3,3,7)
zmap = my_results.zmap_c3;
mask_in = my_results.zmapthresh_clc_c3;
my_title = 'Difference (LEPs-CHEPs)';
plot_zmaps_with_outline(zmap,logical(mask_in),x,y,clims_in,my_title)
%Cz
ax(8) =subplot(3,3,8)
zmap = my_results.zmap_cz;
mask_in = my_results.zmapthresh_clc_cz;
my_title = 'Difference (LEPs-CHEPs)';
plot_zmaps_with_outline(zmap,logical(mask_in),x,y,clims_in,my_title)
%C4
ax(9) =subplot(3,3,9)
zmap = my_results.zmap_c4;
mask_in = my_results.zmapthresh_clc_c4;
my_title = 'Difference (LEPs-CHEPs)';
plot_zmaps_with_outline(zmap,logical(mask_in),x,y,clims_in,my_title)
%colormap(new_color_map_2)
colormap(ax(1),jet);
colormap(ax(2),jet);
colormap(ax(3),jet);
colormap(ax(4),jet);
colormap(ax(5),jet);
colormap(ax(6),jet);
colormap(ax(7),my_custom_map);
colormap(ax(8),my_custom_map);
colormap(ax(9),my_custom_map);



set(gcf,'color','w');



%exportgraphics(f,'CHEPSvsLEPS_final.png','Resolution',600)


%% Testing

vec = [      100;       83;       68;       44;       30;       15;        0];
hex = ['#f8fe34';'#f4771e';'#d94015';'#148b3e';'#085a3f';'#41332d';'#010101'];
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
N = 128;
%N = size(get(gcf,'colormap'),1) % size of the current colormap
map = interp1(vec,raw,linspace(100,0,N),'pchip');

