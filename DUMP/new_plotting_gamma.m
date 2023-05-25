%% New stats script

% load data
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data\round2
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\scripts_functions\Colormaps
load('my_data.mat');
results = load('result_masks.mat');

cd Z:\30_Oscar_Ortiz\LEPs_CHEPs\scripts_functions

%% Getting indeces 
cond1 = 'c1';
cond2 = 'c3';
struct_in1 = cond_1;
struct_in2 = cond_3;
chan_name1 = "C3";
chan_name2 = "C4";
chan_name3 = "Cz";
for k = 1:length(struct_in1.chan_names_stored)
    C_chan = string(struct_in1.chan_names_stored(k));
    if C_chan == lower(chan_name1)
        chan1_IDX=k;
    end
    if C_chan == lower(chan_name2)
        chan2_IDX=k;
    end
    if C_chan == lower(chan_name3)
        chan3_IDX=k;
    end
end

%% clearing s08 of condition 3  since cond 1 is missing 8 (only run 1nce)
temp1 = cond_3.gamma_data([1:7 9:end],:,:,:);
new_mean = squeeze(mean(temp1,1));
cond_3.gamma_data = temp1;
cond_3.gamma_mean = new_mean;

temp1 = cond_3.lf_data([1:7 9:end],:,:,:);
new_mean = squeeze(mean(temp1,1));
cond_3.lf_data = temp1;
cond_3.lf_mean = new_mean;

%%
%Channel 1 - C3
time_window = [-250 800];
v_time = results.g_tfv_time;
plotting_input = 1;
data_in1 = squeeze(cond_3.gamma_data(:,chan1_IDX,:,:));
data_in2 = squeeze(cond_1.gamma_data(:,chan1_IDX,:,:));
mask_sig1 = results.g_c3_zmapthresh_clc_c3;
mask_sig2 = results.g_c1_zmapthresh_clc_c3;
mask_diff = results.gamma_zmapthresh_clc_c3;
[means_output,h,p,ci,stats] = sts_pair_overlap_ttest(data_in1, data_in2,mask_sig1, mask_sig2,mask_diff, time_window, v_time, plotting_input)

%% Channel 2 - C4

data_in1 = squeeze(cond_3.gamma_data(:,chan2_IDX,:,:));
data_in2 = squeeze(cond_1.gamma_data(:,chan2_IDX,:,:));
mask_sig1 = results.g_c3_zmapthresh_clc_c4;
mask_sig2 = results.g_c1_zmapthresh_clc_c4;
mask_diff = results.gamma_zmapthresh_clc_c4;
[means_output,h,p,ci,stats] = sts_pair_overlap_ttest(data_in1, data_in2,mask_sig1, mask_sig2,mask_diff, time_window, v_time, plotting_input)

%Channel 3 - Cz
data_in1 = squeeze(cond_3.gamma_data(:,chan3_IDX,:,:));
data_in2 = squeeze(cond_1.gamma_data(:,chan3_IDX,:,:));
mask_sig1 = results.g_c3_zmapthresh_clc_cz;
mask_sig2 = results.g_c1_zmapthresh_clc_cz;
mask_diff = results.gamma_zmapthresh_clc_cz;
[means_output,h,p,ci,stats] = sts_pair_overlap_ttest(data_in1, data_in2,mask_sig1, mask_sig2,mask_diff, time_window, v_time, plotting_input)




%% Plotting channel 1 - C3 
my_custom_map = fake_parula(250); %inferno fake parula plasma viridis 

figure 
tiledlayout(3,3);
ax(1) = nexttile;
data_in = squeeze(cond_3.gamma_mean(chan1_IDX,:,:));
mask_in =  results.g_c3_zmapthresh_clc_c3;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma C3, LEPS';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);

ax(2) = nexttile;
data_in = squeeze(cond_3.gamma_mean(chan3_IDX,:,:));
mask_in =  results.g_c3_zmapthresh_clc_cz;
x = cond_3.times;
x1 = results.tfv_time_cz;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma Cz, LEPS';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);

ax(3) = nexttile;
data_in = squeeze(cond_3.gamma_mean(chan2_IDX,:,:));
mask_in =  results.g_c3_zmapthresh_clc_c4;
x = cond_3.times;
x1 = results.tfv_time_c4;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma C4, LEPS';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);


% condition 1
ax(4) = nexttile;
data_in = squeeze(cond_1.gamma_mean(chan1_IDX,:,:));
mask_in =  results.g_c1_zmapthresh_clc_c3;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma C3, CHEPs';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);

ax(5) = nexttile;
data_in = squeeze(cond_1.gamma_mean(chan3_IDX,:,:));
mask_in =  results.g_c1_zmapthresh_clc_cz;
x = cond_3.times;
x1 = results.tfv_time_cz;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma Cz, CHEPs';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in)

ax(6)= nexttile;
data_in = squeeze(cond_1.gamma_mean(chan2_IDX,:,:));
mask_in =  results.g_c1_zmapthresh_clc_c4;
x = cond_3.times;
x1 = results.tfv_time_c4;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma C4, CHEPs';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in)

% plot differences
my_clims = [-2.5 2.5];
my_xlims = [-150 800];
%C3
ax(7) = nexttile;
zmap = results.gamma_zmap_c3;
mask_in = results.gamma_zmapthresh_clc_c3;
my_title = 'Difference (LEPs-CHEPs)';
plot_zmaps_with_outline(zmap,logical(mask_in),x1,y,my_clims,my_title)

%Cz
ax(8) = nexttile;
zmap = results.gamma_zmap_cz;
mask_in = results.gamma_zmapthresh_clc_cz;
my_title = 'Difference (LEPs-CHEPs)';
plot_zmaps_with_outline(zmap,logical(mask_in),x1,y,my_clims,my_title)

%Cz
ax(9) = nexttile;
zmap = results.gamma_zmap_c4;
mask_in = results.gamma_zmapthresh_clc_c4;
my_title = 'Difference (LEPs-CHEPs)';
plot_zmaps_with_outline(zmap,logical(mask_in),x1,y,my_clims,my_title)


colormap(ax(1),jet);
colormap(ax(2),jet);
colormap(ax(3),jet);
colormap(ax(4),jet);
colormap(ax(5),jet);
colormap(ax(6),jet);
colormap(ax(7),my_custom_map);
colormap(ax(8),my_custom_map);
colormap(ax(9),my_custom_map);



%%
figure
subplot(3,1,1)
time_window = [-250 800];
v_time = results.g_tfv_time;
plotting_input = 1;
data_in1 = squeeze(cond_3.gamma_data(:,chan1_IDX,:,:));
data_in2 = squeeze(cond_1.gamma_data(:,chan1_IDX,:,:));
mask_sig1 = results.g_c3_zmapthresh_clc_c3;
mask_sig2 = results.g_c1_zmapthresh_clc_c3;
mask_diff = results.gamma_zmapthresh_clc_c3;
nexttile
sts_pair_overlap_ttest(data_in1, data_in2,mask_sig1, mask_sig2,mask_diff, time_window, v_time, plotting_input)


%Channel 3 - Cz
nexttile
data_in1 = squeeze(cond_3.gamma_data(:,chan3_IDX,:,:));
data_in2 = squeeze(cond_1.gamma_data(:,chan3_IDX,:,:));
mask_sig1 = results.g_c3_zmapthresh_clc_cz;
mask_sig2 = results.g_c1_zmapthresh_clc_cz;
mask_diff = results.gamma_zmapthresh_clc_cz;
[means_output,h,p,ci,stats] = sts_pair_overlap_ttest(data_in1, data_in2,mask_sig1, mask_sig2,mask_diff, time_window, v_time, plotting_input);




%% Channel 2 - C4
nexttile
data_in1 = squeeze(cond_3.gamma_data(:,chan2_IDX,:,:));
data_in2 = squeeze(cond_1.gamma_data(:,chan2_IDX,:,:));
mask_sig1 = results.g_c3_zmapthresh_clc_c4;
mask_sig2 = results.g_c1_zmapthresh_clc_c4;
mask_diff = results.gamma_zmapthresh_clc_c4;
[means_output,h,p,ci,stats] = sts_pair_overlap_ttest(data_in1, data_in2,mask_sig1, mask_sig2,mask_diff, time_window, v_time, plotting_input);


%% Trying to concatenate lf as well


all_data = concatenate_all_frequencies(gamma_in,lf_in);
my_results = concatenate_all_results(gamma_results,lf_results);


