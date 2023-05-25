%% Example for how to use npp_stb

load('grouped_data_all_chans_extracted.mat'); % this is the dataset oc the cheps vs leps study. In here we have two structures, one per condition containing the data. 

%% defining parameter
struct_in = cond_1_hand_corrected; % temporary structure to get some channel indeces

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

% some of the input parameters for the function
input_mat = squeeze(struct_in.gamma_data(:,chan1_IDX,:,:)); % input matrix containing all the subjects
n_permutes = 5000; % number of permutations
v_time = struct_in.times; % vector of time
v_freq = struct_in.lf_frex; % vector of frequencies in 
time_window = [-150 800]; % time window that you want to analyze 
plotting_input = 1; % make 1 if you want to plot the responses or zero if you don't. 


%% Using the function


%% condition 1 - gamma
struct_in = cond_1_hand_corrected; % make struct_in be the condition you want to use


% C3 
my_title = ['Cheps low frequencies ' char(chan_name1)];
input_mat = squeeze(struct_in.lf_data(:,chan1_IDX,:,:));
[c1_zmap_c3, c1_zmapthresh_c3, c1_zmapthresh_plc_c3, c1_zmapthresh_clc_c3, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Cheps low frequencies ' char(chan_name2)];
input_mat = squeeze(struct_in.lf_data(:,chan2_IDX,:,:));
[c1_zmap_c4, c1_zmapthresh_c4, c1_zmapthresh_plc_c4, c1_zmapthresh_clc_c4, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Cheps low frequencies ' char(chan_name3)];
input_mat = squeeze(struct_in.lf_data(:,chan3_IDX,:,:));
[c1_zmap_cz, c1_zmapthresh_cz, c1_zmapthresh_plc_cz, c1_zmapthresh_clc_cz, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


%% condition 3  - gamma

struct_in = cond_3_hand_corrected; % make struct_in be the condition you want to use

% C3
my_title = ['Leps low frequencies ' char(chan_name1)];
input_mat = squeeze(struct_in.lf_data(:,chan1_IDX,:,:));
[c3_zmap_c3, c3_zmapthresh_c3, c3_zmapthresh_plc_c3, c3_zmapthresh_clc_c3, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Leps low frequencies ' char(chan_name2)];
input_mat = squeeze(struct_in.lf_data(:,chan2_IDX,:,:));
[c3_zmap_c4, c3_zmapthresh_c4, c3_zmapthresh_plc_c4, c3_zmapthresh_clc_c4, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Leps low frequencies ' char(chan_name3)];
input_mat = squeeze(struct_in.lf_data(:,chan3_IDX,:,:));
[c3_zmap_cz, c3_zmapthresh_cz, c3_zmapthresh_plc_cz, c3_zmapthresh_clc_cz, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


%% Plot cluster corrected responses





