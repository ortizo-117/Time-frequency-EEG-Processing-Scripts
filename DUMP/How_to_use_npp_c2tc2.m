%% do the point by point comparison between cond 1 and 3

addpath Z:/30_Oscar_Ortiz/LEPs_CHEPs/4_grouped_data
load('grouped_data_all_chans_extracted.mat')
struct_in1 = cond_1_clean;
struct_in2 = cond_2_clean;
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


input_mat1 = squeeze(struct_in1.gamma_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in2.gamma_data([1:7 9:end],chan1_IDX,:,:));
n_permutes = 1000;
v_time = struct_in1.times;
v_freq = struct_in1.gamma_frex;
time_window = [-150 800];
plotting_input = 1;

% gamma
input_mat1 = squeeze(struct_in1.gamma_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in2.gamma_data([1:7 9:end],chan1_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 1- C3 gamma';
[gamma_zmap_c3, gamma_zmapthresh_c3, gamma_zmapthresh_plc_c3, gamma_zmapthresh_clc_c3, tfv_time_c3, v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

input_mat1 = squeeze(struct_in1.gamma_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in2.gamma_data([1:7 9:end],chan2_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 1- C4 gamma';
[gamma_zmap_c4, gamma_zmapthresh_c4, gamma_zmapthresh_plc_c4, gamma_zmapthresh_clc_c4, tfv_time_c4, v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



input_mat1 = squeeze(struct_in1.gamma_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in2.gamma_data([1:7 9:end],chan3_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 1- Cz gamma';
[gamma_zmap_cz, gamma_zmapthresh_cz, gamma_zmapthresh_plc_cz, gamma_zmapthresh_clc_cz, tfv_time_cz, v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);
