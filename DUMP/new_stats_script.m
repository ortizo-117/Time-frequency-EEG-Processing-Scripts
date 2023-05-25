%% New stats script
%% Cleaning subjects that are missing conditions...
% condition 1 is missing S08 so the index to remove from condition 3 is 8


cd Z:\30_Oscar_Ortiz\LEPs_CHEPs\scripts_functions
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data\round2

% Loading Data
load ('my_data.mat');

%do for gamma
clean_dat = cond_3.gamma_data([1:7, 9:end],:,:,:); % delete s08
clean_mean = squeeze(mean(clean_dat)); % recompute the mean
cond_3.gamma_data = clean_dat;
cond_3.gamma_mean = clean_mean;

% do for lf
clean_dat = cond_3.lf_data([1:7, 9:end],:,:,:); % delete s08
clean_mean = squeeze(mean(clean_dat)); % recompute the mean
cond_3.lf_data = clean_dat;
cond_3.lf_mean = clean_mean;


struct_in1 = cond_1;
struct_in2 = cond_3;



%% Extraction of low frequency data
struct_in = cond_1;
%% Find indices
chan_name1 = "C3";
chan_name2 = "C4";
chan_name3 = "Cz";
chan_name4 = "Fz";
chan_name5 = "Pz";
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
    if C_chan == lower(chan_name4)
        chan4_IDX=k;
    end
    if C_chan == lower(chan_name5)
        chan5_IDX=k;
    end
end

input_mat = squeeze(struct_in.lf_data(:,chan1_IDX,:,:));
n_permutes = 1000;
v_time = struct_in.times;
v_freq = struct_in.lf_frex;
time_window = [-300 800];
plotting_input = 1;


%% condition 1
% C3
my_title = ['Cheps LF ' char(chan_name1)];
input_mat = squeeze(struct_in.lf_data(:,chan1_IDX,:,:));
[lf_c1_zmap_c3, lf_c1_zmapthresh_c3, lf_c1_zmapthresh_plc_c3, lf_c1_zmapthresh_clc_c3, lf_tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Cheps gamma ' char(chan_name2)];
input_mat = squeeze(struct_in.lf_data(:,chan2_IDX,:,:));
[lf_c1_zmap_c4, lf_c1_zmapthresh_c4, lf_c1_zmapthresh_plc_c4, lf_c1_zmapthresh_clc_c4, lf_tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Cheps gamma ' char(chan_name3)];
input_mat = squeeze(struct_in.lf_data(:,chan3_IDX,:,:));
[lf_c1_zmap_cz, lf_c1_zmapthresh_cz, lf_c1_zmapthresh_plc_cz, lf_c1_zmapthresh_clc_cz, lf_tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



% gamma of fz
struct_in2 = cond_3;
struct_in1 = cond_1;
n_permutes = 2000;
v_time = struct_in1.times;
v_freq = struct_in1.gamma_frex;
time_window = [-300 800];
plotting_input = 1;


my_title =  ['Cheps gamma ' char(chan_name4)];
input_mat = squeeze(struct_in1.gamma_data(:,chan4_IDX,:,:));
[g_c1_zmap_fz, g_c1_zmapthresh_fz, g_c1_zmapthresh_plc_fz, g_c1_zmapthresh_clc_fz, gamma_tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

my_title =  ['Cheps gamma ' char(chan_name5)];
input_mat = squeeze(struct_in1.gamma_data(:,chan5_IDX,:,:));
[g_c1_zmap_pz, g_c1_zmapthresh_pz, g_c1_zmapthresh_plc_pz, g_c1_zmapthresh_clc_pz, gamma_tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



my_title =  ['Leps gamma ' char(chan_name4)];
input_mat = squeeze(struct_in2.gamma_data(:,chan4_IDX,:,:));
[g_c3_zmap_fz, g_c3_zmapthresh_fz, g_c3_zmapthresh_plc_fz, g_c3_zmapthresh_clc_fz, gamma_tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

my_title =  ['Leps gamma ' char(chan_name5)];
input_mat = squeeze(struct_in2.gamma_data(:,chan5_IDX,:,:));
[g_c3_zmap_pz, g_c3_zmapthresh_pz, g_c3_zmapthresh_plc_pz, g_c3_zmapthresh_clc_pz, gamma_tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



%% condition 3
struct_in = cond_3;

% C3
my_title = ['Leps gamma ' char(chan_name1)];
input_mat = squeeze(struct_in.lf_data(:,chan1_IDX,:,:));
[lf_c3_zmap_c3, lf_c3_zmapthresh_c3, lf_c3_zmapthresh_plc_c3, lf_c3_zmapthresh_clc_c3, lf_tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Leps gamma ' char(chan_name2)];
input_mat = squeeze(struct_in.lf_data(:,chan2_IDX,:,:));
[lf_c3_zmap_c4, lf_c3_zmapthresh_c4, lf_c3_zmapthresh_plc_c4, lf_c3_zmapthresh_clc_c4, lf_tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Leps gamma ' char(chan_name3)];
input_mat = squeeze(struct_in.lf_data(:,chan3_IDX,:,:));
[lf_c3_zmap_cz, lf_c3_zmapthresh_cz, lf_c3_zmapthresh_plc_cz, lf_c3_zmapthresh_clc_cz, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%% Running the comparisons 
% load data
results = load('results_updated_V3.mat');
load('my_data.mat')

%%% LF
%C3
struct_in2 = cond_3;
struct_in1 = cond_1;
input_mat1 = squeeze(struct_in1.lf_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in2.lf_data(:,chan1_IDX,:,:));
n_permutes = 1000;
v_time = struct_in1.times;
v_freq = struct_in1.lf_frex;
time_window = [-300 800];
plotting_input = 1;

%%%Condition 1 vs Condition 3
my_title = 'NPP t-test condition 3 - codition 1- C3 lf';
[lowfreq_zmap_c3, lowfreq_zmapthresh_c3, lowfreq_zmapthresh_plc_c3, lowfreq_zmapthresh_clc_c3, lowfreq_tfv_time_c3, lowfreq_v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% Cz
input_mat1 = squeeze(struct_in1.lf_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in2.lf_data(:,chan2_IDX,:,:));
%%%Condition 1 vs Condition 3
my_title = 'NPP t-test condition 3 - codition 1- C4 lf';
[lowfreq_zmap_c4, lowfreq_zmapthresh_c4, lowfreq_zmapthresh_plc_c4, lowfreq_zmapthresh_clc_c4, lowfreq_tfv_time_c4, lowfreq_v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% Cz
input_mat1 = squeeze(struct_in1.lf_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in2.lf_data(:,chan3_IDX,:,:));
%%%Condition 1 vs Condition 3
my_title = 'NPP t-test condition 3 - codition 1- Cz lf';
[lowfreq_zmap_cz, lowfreq_zmapthresh_cz, lowfreq_zmapthresh_plc_cz, lowfreq_zmapthresh_clc_cz, lowfreq_tfv_time_cz, lowfreq_v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


%% Load data
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data\round2
load('my_data.mat')



%% Find indices
chan_name1 = "C3";
chan_name2 = "C4";
chan_name3 = "Cz";
chan_name4 = "Fz";
chan_name5 = "Pz";
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
    if C_chan == lower(chan_name4)
        chan4_IDX=k;
    end
    if C_chan == lower(chan_name5)
        chan5_IDX=k;
    end
end

%% Define input matrices for the statistcs
% C3 vs C4 in the LEPS - gamma
input_mat1 = squeeze(cond_3.gamma_data(:,chan2_IDX,:,:)); % C4
input_mat2 = squeeze(cond_3.gamma_data(:,chan1_IDX,:,:)); % C3
n_permutes = 1000;
v_time = struct_in1.times;
v_freq = struct_in1.gamma_frex;
time_window = [-300 800];
plotting_input = 1;
my_title = 'NPP t-test LEPS c4 vs c3';
[leps_gamma_zmapthresh_c4_vs_c3, leps_gamma_zmapthresh_plc_c4_vs_c3, leps_gamma_zmapthresh_clc_c4_vs_c3, leps_tfv_time_c4_vs_c3, leps_ , leps_v_freq_c4_vs_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


% C3 vs C4 in the CHEPS
input_mat1 = squeeze(cond_1.gamma_data(:,chan2_IDX,:,:)); % C4
input_mat2 = squeeze(cond_1.gamma_data(:,chan1_IDX,:,:)); % C3
n_permutes = 1000;
v_time = struct_in1.times;
v_freq = struct_in1.gamma_frex;
time_window = [-300 800];
plotting_input = 1;
my_title = 'NPP t-test CHEPs c4 vs c3';
[cheps_gamma_zmapthresh_c4_vs_c3, cheps_gamma_zmapthresh_plc_c4_vs_c3, cheps_gamma_zmapthresh_clc_c4_vs_c3, cheps_tfv_time_c4_vs_c3,cheps_ cheps_v_freq_c4_vs_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



%% Running the comparisons 
% load data
results = load('results_updated_V3.mat');
load('my_data.mat')

%%% LF
%C3
struct_in2 = cond_3;
struct_in1 = cond_1;
input_mat1 = squeeze(struct_in1.lf_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in2.lf_data(:,chan1_IDX,:,:));
n_permutes = 1000;
v_time = struct_in1.times;
v_freq = struct_in1.lf_frex;
time_window = [-300 800];
plotting_input = 1;

%%%Condition 1 vs Condition 3

my_title = 'NPP t-test condition 3 - codition 1- C3 lf';
[lowfreq_zmap_c3, lowfreq_zmapthresh_c3, lowfreq_zmapthresh_plc_c3, lowfreq_zmapthresh_clc_c3, lowfreq_tfv_time_c3, lowfreq_v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% Cz
input_mat1 = squeeze(struct_in1.lf_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in2.lf_data(:,chan2_IDX,:,:));
%%%Condition 1 vs Condition 3
my_title = 'NPP t-test condition 3 - codition 1- C4 lf';
[lowfreq_zmap_c4, lowfreq_zmapthresh_c4, lowfreq_zmapthresh_plc_c4, lowfreq_zmapthresh_clc_c4, lowfreq_tfv_time_c4, lowfreq_v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% Cz
input_mat1 = squeeze(struct_in1.lf_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in2.lf_data(:,chan3_IDX,:,:));

%%%Condition 1 vs Condition 3

my_title = 'NPP t-test condition 3 - codition 1- Cz lf';
[lowfreq_zmap_cz, lowfreq_zmapthresh_cz, lowfreq_zmapthresh_plc_cz, lowfreq_zmapthresh_clc_cz, lowfreq_tfv_time_cz, lowfreq_v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);




%% creating overlapping segments
results = load('results_updated_V3.mat');
%load('my_data.mat')

%% plotting
figure 
tiledlayout(3,2);
ax(1) = nexttile; % C3 leps
data_in = squeeze(cond_3.gamma_mean(chan1_IDX,:,:));
mask_in =  results.g_c3_zmapthresh_clc_c3;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma C3, LEPS';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);



ax(2) = nexttile; % C3 leps
data_in = squeeze(cond_1.gamma_mean(chan1_IDX,:,:));
mask_in =  results.g_c1_zmapthresh_clc_c3;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma C4, Cheps';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);

ax(3) = nexttile; % C4 leps
data_in = squeeze(cond_3.gamma_mean(chan2_IDX,:,:));
mask_in =  results.g_c3_zmapthresh_clc_c4;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_1.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma C3, Cheps';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);


ax(4) = nexttile; % C3 cheps
data_in = squeeze(cond_1.gamma_mean(chan2_IDX,:,:));
mask_in =  results.g_c1_zmapthresh_clc_c4;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_1.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma C4, Cheps';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);


% plot differences
my_clims = [-4 4];
my_xlims = [-150 800];
%LEPS
ax(5) = nexttile;
zmap = results.leps_gamma_zmapthresh_c4_vs_c3;
mask_in = results.leps_gamma_zmapthresh_clc_c4_vs_c3;
my_title = 'Difference (C4-C3)';
plot_zmaps_with_outline(zmap,logical(mask_in),x1,y,my_clims,my_title)


%CHEPS
ax(5) = nexttile;
zmap = results.cheps_gamma_zmapthresh_c4_vs_c3;
mask_in = results.cheps_gamma_zmapthresh_clc_c4_vs_c3;
my_title = 'Difference (C4-C3)';
plot_zmaps_with_outline(zmap,logical(mask_in),x1,y,my_clims,my_title)



%% Plotting finalized 
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\scripts_functions\Colormaps

%load('my_data.mat')
results = load('results_updated_V4.mat');

% Join lf and gamma frequencies
cond_1_all = join_frequencies(cond_1);
cond_3_all = join_frequencies(cond_3);
[cond_1_all, cond_3_all]=join_frequencies_2_structures(cond_1,cond_3);
% Do the same for the results
stb_gamma = load('results_stb_gamma.mat');
stb_lf = load('results_stb_lf.mat');
c1vc2_gamma = load('results_c1vc2_gamma.mat');
c1vc2_lf = load('results_c1vc2_lf.mat');


%% Actually plotting
%my_custom_map = fake_parula(250); %inferno fake parula plasma viridis 
load('my_custom_map.mat')
figure

%%%% actually plotting - C3
% Creating masks from results
C3_c1_mask = logical([results.lf_c1_zmapthresh_clc_c3(1:end-1,:);results.g_c1_zmapthresh_clc_c3]);
C3_c3_mask = logical([results.lf_c3_zmapthresh_clc_c3(1:end-1,:);results.g_c3_zmapthresh_clc_c3]);
C3_mask_diff = logical([results.lowfreq_zmapthresh_clc_c3(1:end-1,:);results.gamma_zmapthresh_clc_c3]);
C3_zmap = [results.lowfreq_zmap_c3(1:end-1,:);-results.gamma_zmap_c3];
y = [results.lowfreq_v_freq_c3(1:end-1), results.leps_v_freq_c4_vs_c3];
x = results.tfv_time;

% Cheps
ax(1) = subplot(3,3,1);
data_in = squeeze(cond_1_all.allmean_n(chan1_IDX,:,:));
mask_in =  C3_c1_mask;
x1 = cond_3_all.times;
clims_in = [-1 1];
my_xlims = [min(x) max(x)];
title_in = 'C3, CHEPS';
plot_og_dat_with_mask(data_in,mask_in,x1,x,y,clims_in,my_xlims,title_in);
yline(30,'LineWidth',3,'LineStyle','--','Color','k');
%Leps
ax(2) = subplot(3,3,4);
data_in = squeeze(cond_3_all.allmean_n(chan1_IDX,:,:));
mask_in =  C3_c3_mask;
x1 = cond_3_all.times;
clims_in = [-1 1];
my_xlims = [min(x) max(x)];
title_in = 'C3, LEPS';
plot_og_dat_with_mask(data_in,mask_in,x1,x,y,clims_in,my_xlims,title_in);
yline(30,'LineWidth',3,'LineStyle','--','Color','k');
% Difference
ax(3) = subplot(3,3,7); 
clims_in = [-5 5];
title_in = ["Difference"];
zmap = C3_zmap;
mask_in  = C3_mask_diff;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
yline(30,'LineWidth',3,'LineStyle','--','Color','k');

%%%% actually plotting - Cz
% Creating masks from results
Cz_c1_mask = logical([results.lf_c1_zmapthresh_clc_cz(1:end-1,:);results.g_c1_zmapthresh_clc_cz]);
Cz_c3_mask = logical([results.lf_c3_zmapthresh_clc_cz(1:end-1,:);results.g_c3_zmapthresh_clc_cz]);
Cz_mask_diff = logical([results.lowfreq_zmapthresh_plc_cz(1:end-1,:);results.gamma_zmapthresh_clc_cz]);
Cz_zmap = [results.lowfreq_zmap_cz(1:end-1,:);-results.gamma_zmap_cz];


% Cheps
ax(4) = subplot(3,3,2);
data_in = squeeze(cond_1_all.allmean_n(chan3_IDX,:,:));
mask_in =  Cz_c1_mask;
x1 = cond_3_all.times;
clims_in = [-1 1];
my_xlims = [min(x) max(x)];
title_in = 'Cz, CHEPS';
plot_og_dat_with_mask(data_in,mask_in,x1,x,y,clims_in,my_xlims,title_in);
yline(30,'LineWidth',3,'LineStyle','--','Color','k');
%Leps
ax(5) = subplot(3,3,5);
data_in = squeeze(cond_3_all.allmean_n(chan3_IDX,:,:));
mask_in =  Cz_c3_mask;
x1 = cond_3_all.times;
clims_in = [-1 1];
my_xlims = [min(x) max(x)];
title_in = 'Cz, LEPS';
plot_og_dat_with_mask(data_in,mask_in,x1,x,y,clims_in,my_xlims,title_in);
yline(30,'LineWidth',3,'LineStyle','--','Color','k');

% Difference
ax(6) = subplot(3,3,8);
clims_in = [-5 5];
title_in = ["Difference"];
zmap = Cz_zmap;
mask_in  = Cz_mask_diff;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
yline(30,'LineWidth',3,'LineStyle','--','Color','k');



%%%% actually plotting - Cz
% Creating masks from results
C4_c1_mask = logical([results.lf_c1_zmapthresh_clc_c4(1:end-1,:);results.g_c1_zmapthresh_clc_c4]);
C4_c3_mask = logical([results.lf_c3_zmapthresh_clc_c4(1:end-1,:);results.g_c3_zmapthresh_clc_c4]);
C4_mask_diff = logical([results.lowfreq_zmapthresh_clc_c4(1:end-1,:);results.gamma_zmapthresh_clc_c4]);
C4_zmap = [results.lowfreq_zmap_c4(1:end-1,:);-results.gamma_zmap_c4];


% Cheps
ax(7) = subplot(3,3,3);

data_in = squeeze(cond_1_all.allmean_n(chan2_IDX,:,:));
mask_in =  C4_c1_mask;
x1 = cond_3_all.times;
clims_in = [-1 1];
my_xlims = [min(x) max(x)];
title_in = 'C4, CHEPS';
plot_og_dat_with_mask(data_in,mask_in,x1,x,y,clims_in,my_xlims,title_in);
yline(30,'LineWidth',3,'LineStyle','--','Color','k');

%Leps
ax(8) = subplot(3,3,6);
data_in = squeeze(cond_3_all.allmean_n(chan2_IDX,:,:));
mask_in =  C4_c3_mask;
x1 = cond_3_all.times;
clims_in = [-1 1];
my_xlims = [min(x) max(x)];
title_in = 'Cz, LEPS';
plot_og_dat_with_mask(data_in,mask_in,x1,x,y,clims_in,my_xlims,title_in);
yline(30,'LineWidth',3,'LineStyle','--','Color','k');

% Difference
ax(9) = subplot(3,3,9);

clims_in = [-5 5];
title_in = ["Difference"];
zmap = C4_zmap;
mask_in  = C4_mask_diff;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
yline(30,'LineWidth',3,'LineStyle','--','Color','k');

colormap(ax(1),jet);
colormap(ax(2),jet);
colormap(ax(3),my_custom_map);
colormap(ax(4),jet);
colormap(ax(5),jet);
colormap(ax(6),my_custom_map);
colormap(ax(7),jet);
colormap(ax(8),jet);
colormap(ax(9),my_custom_map);



set(gcf,'color','w');

print -painters -depsc test3.eps

%% pairwise comparisons 
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
[c3_results,c3h,c3p,c3ci,c3stats] = sts_pair_overlap_ttest(data_in1, data_in2,mask_sig1, mask_sig2,mask_diff, time_window, v_time, plotting_input)


%Channel 3 - Cz
nexttile
data_in1 = squeeze(cond_3.gamma_data(:,chan3_IDX,:,:));
data_in2 = squeeze(cond_1.gamma_data(:,chan3_IDX,:,:));
mask_sig1 = results.g_c3_zmapthresh_clc_cz;
mask_sig2 = results.g_c1_zmapthresh_clc_cz;
mask_diff = results.gamma_zmapthresh_clc_cz;
[cz_means_output,czh,czp,czci,czstats] = sts_pair_overlap_ttest(data_in1, data_in2,mask_sig1, mask_sig2,mask_diff, time_window, v_time, plotting_input);




%% Channel 2 - C4
nexttile
data_in1 = squeeze(cond_3.gamma_data(:,chan2_IDX,:,:));
data_in2 = squeeze(cond_1.gamma_data(:,chan2_IDX,:,:));
mask_sig1 = results.g_c3_zmapthresh_clc_c4;
mask_sig2 = results.g_c1_zmapthresh_clc_c4;
mask_diff = results.gamma_zmapthresh_clc_c4;
[c4_means_output,c4h,c4p,c4ci,c4stats] = sts_pair_overlap_ttest(data_in1, data_in2,mask_sig1, mask_sig2,mask_diff, time_window, v_time, plotting_input);



%% Contralateral vs ipsilateral test - gamma

cd Z:\30_Oscar_Ortiz\LEPs_CHEPs\scripts_functions
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data\round2
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\scripts_functions\Colormaps


% Loading Data
load ('my_data.mat');

%do for gamma
clean_dat = cond_3.gamma_data([1:7, 9:end],:,:,:); % delete s08
clean_mean = squeeze(mean(clean_dat)); % recompute the mean
cond_3.gamma_data = clean_dat;
cond_3.gamma_mean = clean_mean;

% do for lf
clean_dat = cond_3.lf_data([1:7, 9:end],:,:,:); % delete s08
clean_mean = squeeze(mean(clean_dat)); % recompute the mean
cond_3.lf_data = clean_dat;
cond_3.lf_mean = clean_mean;


struct_in1 = cond_1;
struct_in2 = cond_3;


%load('my_data.mat')
results = load('results_updated_V4.mat');


struct_in = cond_1;
chan_name1 = "C3";
chan_name2 = "C4";
chan_name3 = "Cz";
for k = 1:length(struct_in.chan_names_stored)
    C_chan = string(struct_in.chan_names_stored(k));
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




%% Inputs

plotting_input = 1;
time_w = [0 600]; % subregion of time you are interested in
freq_w = [30 100]; % subregion of frequency you are interested in
v_time = cond_1.times; % vector of time matching the size of data
v_freq = cond_1.gamma_frex; % vector of frequencies matching the size of data
v_time_results = results.tfv_time; %subset time for results

data_in1 = squeeze(cond_3.gamma_data(:,chan1_IDX,:,:)); % C3
data_in2 = squeeze(cond_3.gamma_data(:,chan2_IDX,:,:)); % C4
mask_sig1 = logical(results.g_c3_zmapthresh_clc_c3); % significance mask
mask_sig2 = logical(results.g_c3_zmapthresh_clc_c4); % significance mask

title_in = 'LEPs ';
norp = 'negative'; % wether you want to look at negative or positive deflections
[LEPS_ERD_means_output,LEPS_ERD_npix_output,LEPS_ERD_w_means_output,LEPS_ERD_means_results,LEPS_ERD_n_pix_results,LEPS_ERD_w_mean_results] = sts_roi_comparisons(title_in,data_in1, data_in2,mask_sig1, mask_sig2,norp, time_w,freq_w,v_time,v_freq,v_time_results,plotting_input)
norp = 'positive'; % wether you want to look at negative or positive deflections
[LEPS_ERS_means_output,LEPS_ERS_npix_output,LEPS_ERS_w_means_output,LEPS_ERS_means_results,LEPS_ERS_n_pix_results,LEPS_ERS_w_mean_results] = sts_roi_comparisons(title_in,data_in1, data_in2,mask_sig1, mask_sig2,norp, time_w,freq_w,v_time,v_freq,v_time_results,plotting_input)



data_in1 = squeeze(cond_1.gamma_data(:,chan1_IDX,:,:)); % C3
data_in2 = squeeze(cond_1.gamma_data(:,chan2_IDX,:,:)); % C4
mask_sig1 = logical(results.g_c1_zmapthresh_clc_c3); % significance mask
mask_sig2 = logical(results.g_c1_zmapthresh_clc_c4); % significance mask

title_in = 'CHEPs ';
norp = 'negative'; % wether you want to look at negative or positive deflections
[CHEPS_ERD_means_output,CHEPS_ERD_npix_output,CHEPS_ERD_w_means_output,CHEPS_ERD_means_results,CHEPS_ERD_n_pix_results,CHEPS_ERD_w_mean_results] = sts_roi_comparisons(title_in,data_in1, data_in2,mask_sig1, mask_sig2,norp, time_w,freq_w,v_time,v_freq,v_time_results,plotting_input)
norp = 'positive'; % wether you want to look at negative or positive deflections
[CHEPS_ERS_means_output,CHEPS_ERS_npix_output,CHEPS_ERS_w_means_output,CHEPS_ERS_means_results,CHEPS_ERS_n_pix_results,CHEPS_ERS_w_mean_results] = sts_roi_comparisons(title_in,data_in1,data_in2,mask_sig1, mask_sig2,norp, time_w,freq_w,v_time,v_freq,v_time_results,plotting_input)


%% doing by selection of ROI
v_time = cond_1.times; % vector of time matching the size of data
v_freq = cond_1.gamma_frex; % vector of frequencies matching the size of data
v_time_results = results.tfv_time; %subset time for results



data_in1 = squeeze(cond_3.gamma_data(:,chan1_IDX,:,:)); % C3
data_in2 = squeeze(cond_3.gamma_data(:,chan2_IDX,:,:)); % C4
mask_sig1 = logical(results.g_c3_zmapthresh_clc_c3); % significance mask
mask_sig2 = logical(results.g_c3_zmapthresh_clc_c4); % significance mask
title_in = 'LEPs ';
[LEPS_output,LEPS_f1,LEPS_f2] = sts_roi_comparisons_hand_selected(title_in,data_in1, data_in2,mask_sig1, mask_sig2,v_time,v_freq,v_time_results);


data_in1 = squeeze(cond_1.gamma_data(:,chan1_IDX,:,:)); % C3
data_in2 = squeeze(cond_1.gamma_data(:,chan2_IDX,:,:)); % C4
mask_sig1 = logical(results.g_c1_zmapthresh_clc_c3); % significance mask
mask_sig2 = logical(results.g_c1_zmapthresh_clc_c4); % significance mask
title_in = 'CHEPs ';
[CHEPS_output,CHEPS_f1,CHEPS_f2] = sts_roi_comparisons_hand_selected(title_in,data_in1, data_in2,mask_sig1, mask_sig2,v_time,v_freq,v_time_results);


LEPS_ERS_w_mean_results.h 
CHEPS_ERS_w_mean_results.h 
LEPS_ERS_w_mean_results.stats.tstat
CHEPS_ERS_w_mean_results.stats.tstat 
LEPS_ERS_w_mean_results.p 
CHEPS_ERS_w_mean_results.p


LEPS_ERD_w_mean_results.h 
CHEPS_ERD_w_mean_results.h 
LEPS_ERD_w_mean_results.stats.tstat
CHEPS_ERD_w_mean_results.stats.tstat 
LEPS_ERD_w_mean_results.p 
CHEPS_ERD_w_mean_results.p



LEPS_ERS_n_pix_results.h 
CHEPS_ERS_n_pix_results.h 
LEPS_ERS_n_pix_results.stats.tstat
CHEPS_ERS_n_pix_results.stats.tstat 
LEPS_ERS_n_pix_results.p 
CHEPS_ERS_n_pix_results.p


LEPS_ERD_n_pix_results.h 
CHEPS_ERD_n_pix_results.h 
LEPS_ERD_n_pix_results.stats.tstat
CHEPS_ERD_n_pix_results.stats.tstat 
LEPS_ERD_n_pix_results.p 
CHEPS_ERD_n_pix_results.p




LEPS_ERS_means_results.h 
CHEPS_ERS_means_results.h 
LEPS_ERS_means_results.stats.tstat
CHEPS_ERS_means_results.stats.tstat 
LEPS_ERS_means_results.p 
CHEPS_ERS_means_results.p


LEPS_ERD_means_results.h 
CHEPS_ERD_means_results.h 
LEPS_ERD_means_results.stats.tstat
CHEPS_ERD_means_results.stats.tstat 
LEPS_ERD_means_results.p 
CHEPS_ERD_means_results.p

%% Plotting stuff 
results = load('results_V5.mat');
load('my_data.mat');

struct_in1 = cond_1;

chan_name1 = "C3";
chan_name2 = "C4";
chan_name3 = "Cz";
chan_name4 = "Fz";
chan_name5 = "Pz";
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
    if C_chan == lower(chan_name4)
        chan4_IDX=k;
    end
    if C_chan == lower(chan_name5)
        chan5_IDX=k;
    end
end



% cond 1 
% C3
figure
data_in = squeeze(cond_1.gamma_mean(chan1_IDX,:,:));
mask_in =  results.g_c1_zmapthresh_clc_c3;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma C3, CHEPS';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);
set(gcf,'color','w');
exportgraphics(gcf,'cheps_c3.eps','ContentType','vector')



% C4
figure
data_in = squeeze(cond_1.gamma_mean(chan2_IDX,:,:));
mask_in =  results.g_c1_zmapthresh_clc_c4;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma C4, CHEPS';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);
set(gcf,'color','w');
exportgraphics(gcf,'cheps_c4.eps','ContentType','vector')


% Cz
figure
data_in = squeeze(cond_1.gamma_mean(chan3_IDX,:,:));
mask_in =  results.g_c1_zmapthresh_clc_cz;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma Cz, CHEPS';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);
set(gcf,'color','w');
exportgraphics(gcf,'cheps_cz.eps','ContentType','vector')



% Fz
figure
data_in = squeeze(cond_1.gamma_mean(chan4_IDX,:,:));
mask_in =  results.g_c1_zmapthresh_clc_fz;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma Fz, CHEPS';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);
set(gcf,'color','w');
exportgraphics(gcf,'cheps_fz.eps','ContentType','vector')



% Pz
figure
data_in = squeeze(cond_1.gamma_mean(chan5_IDX,:,:));
mask_in =  results.g_c1_zmapthresh_clc_pz;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma Pz, CHEPS';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);
set(gcf,'color','w');
exportgraphics(gcf,'cheps_pz.eps','ContentType','vector')




% cond 3
% C3
figure
data_in = squeeze(cond_3.gamma_mean(chan1_IDX,:,:));
mask_in =  results.g_c3_zmapthresh_clc_c3;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma C3, LEPs';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);
set(gcf,'color','w');
exportgraphics(gcf,'leps_c3.eps','ContentType','vector')


% C4
figure
data_in = squeeze(cond_3.gamma_mean(chan2_IDX,:,:));
mask_in =  results.g_c3_zmapthresh_clc_c4;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma C4, LEPS';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);
set(gcf,'color','w');
exportgraphics(gcf,'leps_c4.eps','ContentType','vector')


% Cz
figure
data_in = squeeze(cond_3.gamma_mean(chan3_IDX,:,:));
mask_in =  results.g_c3_zmapthresh_clc_cz;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma Cz, LEPS';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);
set(gcf,'color','w');
exportgraphics(gcf,'leps_cz.eps','ContentType','vector')



% Fz
figure
data_in = squeeze(cond_3.gamma_mean(chan4_IDX,:,:));
mask_in =  results.g_c3_zmapthresh_clc_fz;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma Fz, LEPS';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);
set(gcf,'color','w');
exportgraphics(gcf,'leps_fz.eps','ContentType','vector')


% Pz
figure
data_in = squeeze(cond_3.gamma_mean(chan5_IDX,:,:));
mask_in =  results.g_c3_zmapthresh_clc_pz;
x = cond_3.times;
x1 = results.tfv_time_c3;
y = cond_3.gamma_frex;
clims_in = [-2 2];
my_xlims = [min(x1) max(x1)];
title_in = 'Gamma Pz, LEPS';
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,clims_in,my_xlims,title_in);
set(gcf,'color','w');
exportgraphics(gcf,'leps_pz.eps','ContentType','vector')

