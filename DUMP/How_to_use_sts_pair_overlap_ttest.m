the 

%% How to use sts_pair_overlap_ttest

% Loading in data
gamma_in = load('nppc_gamma.mat');
my_results  = load('npp_results.mat');
load('grouped_data_all_chans_extracted.mat');


%% getting time parameters
time_window = [-150 800];
v_time = gamma_in.x;

time_s = dsearchn(v_time',time_window(1)); % start time for tf window
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
tfv_time  = v_time(time_s:time_e);
nTimepoints = numel(tfv_time);
chan_index = 1;

% Defining the data
data_in1 = squeeze(cond_1_hand_corrected.gamma_data(:,chan_index,:,time_s:time_e));
data_in2 = squeeze(cond_3_hand_corrected.gamma_data([1:7 9:end],chan_index,:,time_s:time_e)); % this one has weird subsetting cause i had un even number of subjects per condition

% Defining the threshold masks
mask_sig1 = gamma_in.c1_zmapthresh_clc_c3; % mask for significant differences in condition 1
mask_sig2 = gamma_in.c3_zmapthresh_clc_c3; % mask for significant differences in condition 2
mask_diff = my_results.gamma_zmapthresh_clc_c3~=0; % mask for where the two conditions where significantly different

%% Using the function
plotting_input = 1;
[means_output,h,p,ci,stats] = sts_pair_overlap_ttest(data_in1, data_in2,mask_sig1, mask_sig2,mask_diff, time_window, v_time, plotting_input);



