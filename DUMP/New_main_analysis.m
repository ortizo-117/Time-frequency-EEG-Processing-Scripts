%% Indicate path with all folders 
clear all 
close all

% addpath C:\Users\Kip\Documents\eeglab2021.1 % add eeglab path
cd  Z:\30_Oscar_Ortiz\LEPs_CHEPs\scripts_functions
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\Raw
addpath C:\Users\User\Documents\MATLAB\eeglab2021.0 % for oscar at home
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\Raw\Participant_C06_Trial 1,3 & cond 2, 4'; % folder 
addpath C:\Users\Kip\Documents\eeglab2021.1

%% preprocess
saving_folder = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\1_preprocessed_data';
preproces_EEG_folder(pathname, saving_folder,'*epoch_pruned.set')

%% refiltering for sta
saving_folder = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\sta_files\processed_files';
extension_to_look_for = '*low_freq.set';
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\sta_files';
filter_ranges = [0.1 30];
saving_extension = '_lowpass_for_sta';
 filter_files(pathname,saving_folder,extension_to_look_for,saving_extension,filter_ranges)

%% blink and artifact detection step
addpath C:\Users\Kip\Documents\eeglab2021.1;
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\1_preprocessed_data';
saving_folder = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\2_cleaned_data';
%manual_blink_rejection(pathname);

file_struct_list = dir([pathname filesep() '*Participant_C04_Trial_4_Ear_ref_filt_ICA_epoch.set']);  %% get list of .set files in the pathname specified

filename_cell_list = {file_struct_list.name};  %% extract the filenames into a cellarray

filename_list=deblank(char(filename_cell_list));

endout=regexp(pathname,filesep,'split');
saving_name = [char(endout(end-1)) '_' char(endout(end))];

length_filename=size(filename_list);
 
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;  
 
 for k = 1:length_filename(1)
     filename=deblank(filename_list(k, :));   
     [~,filename_text,~]= fileparts(strcat(pathname,'\',filename));   
      EEG.etc.eeglabvers = '14.1.1'; % this tracks which version of EEGLAB is being used, you may ignore it
      EEG = pop_loadset('filename',filename,'filepath',pathname); %loads the specified file into eeglab
      %Inspsection of components by map
      EEG.etc.eeglabvers = '2021.1'; % this tracks which version of EEGLAB is being used, you may ignore it
     
     j=0;
     while j == 0 % in case you mess up a component selection you can try again
     % Inspect components reflecting blinks and noise visually
     pop_selectcomps(EEG, [1:20] );
     pause 
     prompt = 'Are you happy with your selection of components? (input 1 to continue or 0 to do again)';
     j = input(prompt);
     end
    
%      pop_selectcomps(EEG,[1:7] );
%      pause
%       prompt = 'Blink Component selection ok? (1 = yes, continue. 0 = no, do again)';
%       x = input(prompt);
%       if x == 0
%          pop_selectcomps(EEG, [1:7] );
%          x = input(prompt);
%       else
%       end  
      % pause
       EEG = pop_subcomp( EEG,[], 0,1);
       EEG = eeg_checkset( EEG );
       
       m=0;
       fh1 = figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
       while m == 0;
         pop_eegplot( EEG, 0, 1, 1);
         pause 
         prompt = 'Are you happy with your deletion of segments? (input 1 to continue or 0 to do again)';
         m = input(prompt);
       end
       
       try
       close(fh1);
       catch
       end
       
        saving_name = [filename_text '_pruned'];
        EEG = pop_saveset( EEG, 'filename',saving_name,'filepath', saving_folder);
        EEG = eeg_checkset( EEG );
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        eeglab redraw;
        fprintf('\n\n\n %i percent processing folder \n\n\n',k/length_filename(1)*100);      
 end
 
 
%  trial = 'Participant_C07_Trial_3_Ear_ref_filt_ICA_epoch.set'
%  index_i = 0;
%  for i = 1:length(filename_cell_list)
%      if string(trial) == string(filename_cell_list(i))
%          index_i = i;
%          break 
%      end
%  end
%  
 
%% time-frequency decomposition 

pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\2_cleaned_data'; % folder 

files_not_processed = gamma_detection(pathname);


%% Grouping of time-frequency data
addpath Z:\29_Hannah_Goodings\1_Projects\CAP_LONG\CHEP_Only\3_Data_Export
my_chans_in = {'c3','cz','c4','fz','pz'};
% Condition 1
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\gamma_all_extracted\Cond1';
cond_1 = grouping_data_specify_chans(pathname,my_chans_in);

% Condition 2
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\gamma_all_extracted\Cond2';
cond_2 =  grouping_data_specify_chans(pathname,my_chans_in);

% Condition 3
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\gamma_all_extracted\Cond3';
cond_3 =  grouping_data_specify_chans(pathname,my_chans_in);

% Condition 4
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\gamma_all_extracted\Cond4';
cond_4 =  grouping_data_specify_chans(pathname,my_chans_in);


%% Getting only the subjects that are available for all conditions 
% Missing:
%  Cond 1 - Participant 8
%  Cond 2 - Participant 8 and 9 
%  Cond 3 - None Missing
%  Cond 4 - None Missing
% Participant 3 should also be removed for bad data quality
% Overall, we should remove participants 3, 8, 9 
temp = size(cond_1.filenamelist);
cond_length = temp(1);
keep_IDX_Cond1 = [1:2 4:7 9:cond_length]; % need to remove 9

temp = size(cond_2.filenamelist);
cond_length = temp(1);
keep_IDX_Cond2 = [1:2 4:cond_length]; % not removing anything 

temp = size(cond_3.filenamelist);
cond_length = temp(1);
keep_IDX_Cond3 = [1:2 4:7 10:cond_length]; % need to remove 8 and 9 

temp = size(cond_4.filenamelist);
cond_length = temp(1);
keep_IDX_Cond4 = [1:2 4:7 10:cond_length]; % need to remove 8 and 9


cond_1_clean = cond_1;
cond_2_clean = cond_2;
cond_3_clean = cond_3;
cond_4_clean = cond_4;

gamma_data_in = cond_1.gamma_data;
lf_data_in = cond_1.lf_data;
cond_1_clean.gamma_data = gamma_data_in(keep_IDX_Cond1,:,:,:);
cond_1_clean.lf_data = lf_data_in(keep_IDX_Cond1,:,:,:);
cond_1_clean.gamma_mean = squeeze(mean(gamma_data_in(keep_IDX_Cond1,:,:,:)));
cond_1_clean.lf_mean = squeeze(mean(lf_data_in(keep_IDX_Cond1,:,:,:)));


gamma_data_in = cond_2.gamma_data;
lf_data_in = cond_2.lf_data;
cond_2_clean.gamma_data = gamma_data_in(keep_IDX_Cond2,:,:,:);
cond_2_clean.lf_data = lf_data_in(keep_IDX_Cond2,:,:,:);
cond_2_clean.gamma_mean = squeeze(mean(gamma_data_in(keep_IDX_Cond2,:,:,:)));
cond_2_clean.lf_mean = squeeze(mean(lf_data_in(keep_IDX_Cond2,:,:,:)));

gamma_data_in = cond_3.gamma_data;
lf_data_in = cond_3.lf_data;
cond_3_clean.gamma_data = gamma_data_in(keep_IDX_Cond3,:,:,:);
cond_3_clean.lf_data = lf_data_in(keep_IDX_Cond3,:,:,:);
cond_3_clean.gamma_mean = squeeze(mean(gamma_data_in(keep_IDX_Cond3,:,:,:)));
cond_3_clean.lf_mean = squeeze(mean(lf_data_in(keep_IDX_Cond3,:,:,:)));




gamma_data_in = cond_4.gamma_data;
lf_data_in = cond_4.lf_data;
cond_4_clean.gamma_data = gamma_data_in(keep_IDX_Cond4,:,:,:);
cond_4_clean.lf_data = lf_data_in(keep_IDX_Cond4,:,:,:);
cond_4_clean.gamma_mean = squeeze(mean(gamma_data_in(keep_IDX_Cond4,:,:,:)));
cond_4_clean.lf_mean = squeeze(mean(lf_data_in(keep_IDX_Cond4,:,:,:)));

%% Saving new data
saving_folder = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data';
saving_name = 'data_grouped_new.mat';
save([saving_folder filesep() saving_name],'cond_1_clean','cond_2_clean','cond_3_clean','cond_4_clean')


%%  Non parametric permutation using the function
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data
load('data_grouped_new.mat')
struct_in = cond_1_clean;

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

%%Gamma npp
input_mat = squeeze(struct_in.gamma_data(:,chan1_IDX,:,:));
n_permutes = 5000;
v_time = struct_in.times;
v_freq = struct_in.gamma_frex;
time_window = [-150 800];
plotting_input = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% condition 1 Gamma
struct_in = cond_1_clean;
% C3
my_title = ['Baseline Cheps gamma ' char(chan_name1)];
input_mat = squeeze(struct_in.gamma_data(:,chan1_IDX,:,:));
[c1_zmap_c3, c1_zmapthresh_c3, c1_zmapthresh_plc_c3, c1_zmapthresh_clc_c3, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Baseline Cheps gamma ' char(chan_name2)];
input_mat = squeeze(struct_in.gamma_data(:,chan2_IDX,:,:));
[c1_zmap_c4, c1_zmapthresh_c4, c1_zmapthresh_plc_c4, c1_zmapthresh_clc_c4, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Baseline Cheps gamma ' char(chan_name3)];
input_mat = squeeze(struct_in.gamma_data(:,chan3_IDX,:,:));
[c1_zmap_cz, c1_zmapthresh_cz, c1_zmapthresh_plc_cz, c1_zmapthresh_clc_cz, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% condition 3 Gamma
struct_in = cond_3_clean;
% C3
my_title = ['Baseline Leps gamma ' char(chan_name1)];
input_mat = squeeze(struct_in.gamma_data(:,chan1_IDX,:,:));
[c3_zmap_c3, c3_zmapthresh_c3, c3_zmapthresh_plc_c3, c3_zmapthresh_clc_c3, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Baseline Leps gamma ' char(chan_name2)];
input_mat = squeeze(struct_in.gamma_data(:,chan2_IDX,:,:));
[c3_zmap_c4, c3_zmapthresh_c4, c3_zmapthresh_plc_c4, c3_zmapthresh_clc_c4, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Baseline Leps gamma ' char(chan_name3)];
input_mat = squeeze(struct_in.gamma_data(:,chan3_IDX,:,:));
[c3_zmap_cz, c3_zmapthresh_cz, c3_zmapthresh_plc_cz, c3_zmapthresh_clc_cz, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% condition 2 Gamma
struct_in = cond_2_clean;
% C3
my_title = ['Cap Cheps gamma ' char(chan_name1)];
input_mat = squeeze(struct_in.gamma_data(:,chan1_IDX,:,:));
[c2_zmap_c3, c2_zmapthresh_c3, c2_zmapthresh_plc_c3, c2_zmapthresh_clc_c3, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Cap Cheps gamma ' char(chan_name2)];
input_mat = squeeze(struct_in.gamma_data(:,chan2_IDX,:,:));
[c2_zmap_c4, c2_zmapthresh_c4, c2_zmapthresh_plc_c4, c2_zmapthresh_clc_c4, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Cap Cheps gamma ' char(chan_name3)];
input_mat = squeeze(struct_in.gamma_data(:,chan3_IDX,:,:));
[c2_zmap_cz, c2_zmapthresh_cz, c2_zmapthresh_plc_cz, c2_zmapthresh_clc_cz, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% condition 4 Gamma
struct_in = cond_4_clean;
% C3
my_title = ['Cap Leps gamma ' char(chan_name1)];
input_mat = squeeze(struct_in.gamma_data(:,chan1_IDX,:,:));
[c4_zmap_c3, c4_zmapthresh_c3, c4_zmapthresh_plc_c3, c4_zmapthresh_clc_c3, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Cap Leps gamma ' char(chan_name2)];
input_mat = squeeze(struct_in.gamma_data(:,chan2_IDX,:,:));
[c4_zmap_c4, c4_zmapthresh_c4, c4_zmapthresh_plc_c4, c4_zmapthresh_clc_c4, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Cap Leps gamma ' char(chan_name3)];
input_mat = squeeze(struct_in.gamma_data(:,chan3_IDX,:,:));
[c4_zmap_cz, c4_zmapthresh_cz, c4_zmapthresh_plc_cz, c4_zmapthresh_clc_cz, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%% LF nppc
input_mat = squeeze(struct_in.lf_data(:,chan1_IDX,:,:));
n_permutes = 5000;
v_time = struct_in.times;
v_freq = struct_in.lf_frex;
time_window = [-150 800];
plotting_input = 1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% condition 1 lf
struct_in = cond_1_clean;
% C3
my_title = ['Baseline Cheps lf ' char(chan_name1)];
input_mat = squeeze(struct_in.lf_data(:,chan1_IDX,:,:));
[c1_zmap_c3, c1_zmapthresh_c3, c1_zmapthresh_plc_c3, c1_zmapthresh_clc_c3, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Baseline Cheps lf ' char(chan_name2)];
input_mat = squeeze(struct_in.lf_data(:,chan2_IDX,:,:));
[c1_zmap_c4, c1_zmapthresh_c4, c1_zmapthresh_plc_c4, c1_zmapthresh_clc_c4, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Baseline Cheps ;f ' char(chan_name3)];
input_mat = squeeze(struct_in.lf_data(:,chan3_IDX,:,:));
[c1_zmap_cz, c1_zmapthresh_cz, c1_zmapthresh_plc_cz, c1_zmapthresh_clc_cz, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% condition 3 Gamma
struct_in = cond_3_clean;
% C3
my_title = ['Baseline Leps lf ' char(chan_name1)];
input_mat = squeeze(struct_in.lf_data(:,chan1_IDX,:,:));
[c3_zmap_c3, c3_zmapthresh_c3, c3_zmapthresh_plc_c3, c3_zmapthresh_clc_c3, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Baseline Leps lf ' char(chan_name2)];
input_mat = squeeze(struct_in.lf_data(:,chan2_IDX,:,:));
[c3_zmap_c4, c3_zmapthresh_c4, c3_zmapthresh_plc_c4, c3_zmapthresh_clc_c4, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Baseline Leps lf ' char(chan_name3)];
input_mat = squeeze(struct_in.lf_data(:,chan3_IDX,:,:));
[c3_zmap_cz, c3_zmapthresh_cz, c3_zmapthresh_plc_cz, c3_zmapthresh_clc_cz, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% condition 2 Gamma
struct_in = cond_2_clean;
% C3
my_title = ['Cap Cheps lf ' char(chan_name1)];
input_mat = squeeze(struct_in.lf_data(:,chan1_IDX,:,:));
[c2_zmap_c3, c2_zmapthresh_c3, c2_zmapthresh_plc_c3, c2_zmapthresh_clc_c3, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Cap Cheps lf ' char(chan_name2)];
input_mat = squeeze(struct_in.lf_data(:,chan2_IDX,:,:));
[c2_zmap_c4, c2_zmapthresh_c4, c2_zmapthresh_plc_c4, c2_zmapthresh_clc_c4, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Cap Cheps lf ' char(chan_name3)];
input_mat = squeeze(struct_in.lf_data(:,chan3_IDX,:,:));
[c2_zmap_cz, c2_zmapthresh_cz, c2_zmapthresh_plc_cz, c2_zmapthresh_clc_cz, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% condition 4 Gamma
struct_in = cond_4_clean;
% C3
my_title = ['Cap Leps lf ' char(chan_name1)];
input_mat = squeeze(struct_in.lf_data(:,chan1_IDX,:,:));
[c4_zmap_c3, c4_zmapthresh_c3, c4_zmapthresh_plc_c3, c4_zmapthresh_clc_c3, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Cap Leps lf ' char(chan_name2)];
input_mat = squeeze(struct_in.lf_data(:,chan2_IDX,:,:));
[c4_zmap_c4, c4_zmapthresh_c4, c4_zmapthresh_plc_c4, c4_zmapthresh_clc_c4, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Cap Leps lf ' char(chan_name3)];
input_mat = squeeze(struct_in.lf_data(:,chan3_IDX,:,:));
[c4_zmap_cz, c4_zmapthresh_cz, c4_zmapthresh_plc_cz, c4_zmapthresh_clc_cz, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


%% Point by point comparisons

%% do the point by point comparison between cond 1 and 3

addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data
load('data_grouped_new.mat')
struct_in1 = cond_1_clean;
struct_in3 = cond_3_clean;
struct_in2 = cond_2_clean;
struct_in4 = cond_4_clean;


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

%%%%%%%%%%%%%%%%%%% Gamma
%%%Condition 1 vs Condition 3
input_mat1 = squeeze(struct_in1.gamma_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in3.gamma_data(:,chan1_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 1- C3 gamma';
[gamma_zmap_c3, gamma_zmapthresh_c3, gamma_zmapthresh_plc_c3, gamma_zmapthresh_clc_c3, tfv_time_c3, v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

input_mat1 = squeeze(struct_in1.gamma_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in3.gamma_data(:,chan2_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 1- C4 gamma';
[gamma_zmap_c4, gamma_zmapthresh_c4, gamma_zmapthresh_plc_c4, gamma_zmapthresh_clc_c4, tfv_time_c4, v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



input_mat1 = squeeze(struct_in1.gamma_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in3.gamma_data(:,chan3_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 1- C4 gamma';
[gamma_zmap_cz, gamma_zmapthresh_cz, gamma_zmapthresh_plc_cz, gamma_zmapthresh_clc_cz, tfv_time_cz, v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);




% lf
v_freq = struct_in1.lf_frex;
input_mat1 = squeeze(struct_in1.lf_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in3.lf_data(:,chan1_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 1- C3 low frequency';
[lf_zmap_c3, lf_zmapthresh_c3, lf_zmapthresh_plc_c3, lf_zmapthresh_clc_c3, tfv_time_c3, lf_v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

input_mat1 = squeeze(struct_in1.lf_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in3.lf_data(:,chan2_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 1- C4 low frequency';
[lf_zmap_c4, lf_zmapthresh_c4, lf_zmapthresh_plc_c4, lf_zmapthresh_clc_c4, tfv_time_c4, lf_v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


input_mat1 = squeeze(struct_in1.lf_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in3.lf_data(:,chan3_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 1- C4 low frequency';
[lf_zmap_cz, lf_zmapthresh_cz, lf_zmapthresh_plc_cz, lf_zmapthresh_clc_cz, tfv_time_cz, lf_v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



%%%%%%%% Condition 1 vs Condition 2

n_permutes = 1000;
v_time = struct_in1.times;
v_freq = struct_in1.gamma_frex;
time_window = [-150 800];
plotting_input = 1;
input_mat1 = squeeze(struct_in1.gamma_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in2.gamma_data(:,chan1_IDX,:,:));
my_title = 'NPP t-test condition 2 - codition 1- C3 gamma';
[gamma_zmap_c3, gamma_zmapthresh_c3, gamma_zmapthresh_plc_c3, gamma_zmapthresh_clc_c3, tfv_time_c3, v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

input_mat1 = squeeze(struct_in1.gamma_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in2.gamma_data(:,chan2_IDX,:,:));
my_title = 'NPP t-test condition 2 - codition 1- C4 gamma';
[gamma_zmap_c4, gamma_zmapthresh_c4, gamma_zmapthresh_plc_c4, gamma_zmapthresh_clc_c4, tfv_time_c4, v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



input_mat1 = squeeze(struct_in1.gamma_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in2.gamma_data(:,chan3_IDX,:,:));
my_title = 'NPP t-test condition 2 - codition 1- C4 gamma';
[gamma_zmap_cz, gamma_zmapthresh_cz, gamma_zmapthresh_plc_cz, gamma_zmapthresh_clc_cz, tfv_time_cz, v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);




% lf
v_freq = struct_in1.lf_frex;
input_mat1 = squeeze(struct_in1.lf_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in2.lf_data(:,chan1_IDX,:,:));
my_title = 'NPP t-test condition 2 - codition 1- C3 low frequency';
[lf_zmap_c3, lf_zmapthresh_c3, lf_zmapthresh_plc_c3, lf_zmapthresh_clc_c3, tfv_time_c3, lf_v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

input_mat1 = squeeze(struct_in1.lf_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in2.lf_data(:,chan2_IDX,:,:));
my_title = 'NPP t-test condition 2 - codition 1- C4 low frequency';
[lf_zmap_c4, lf_zmapthresh_c4, lf_zmapthresh_plc_c4, lf_zmapthresh_clc_c4, tfv_time_c4, lf_v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


input_mat1 = squeeze(struct_in1.lf_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in2.lf_data(:,chan3_IDX,:,:));
my_title = 'NPP t-test condition 2 - codition 1- C4 low frequency';
[lf_zmap_cz, lf_zmapthresh_cz, lf_zmapthresh_plc_cz, lf_zmapthresh_clc_cz, tfv_time_cz, lf_v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



%%%%%%%%%%%%%%%%% Cond 1 vs Cond 4

n_permutes = 1000;
v_time = struct_in1.times;
v_freq = struct_in1.gamma_frex;
time_window = [-150 800];
plotting_input = 1;
input_mat1 = squeeze(struct_in1.gamma_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in4.gamma_data(:,chan1_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 1- C3 gamma';
[gamma_zmap_c3, gamma_zmapthresh_c3, gamma_zmapthresh_plc_c3, gamma_zmapthresh_clc_c3, tfv_time_c3, v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

input_mat1 = squeeze(struct_in1.gamma_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in4.gamma_data(:,chan2_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 1- C4 gamma';
[gamma_zmap_c4, gamma_zmapthresh_c4, gamma_zmapthresh_plc_c4, gamma_zmapthresh_clc_c4, tfv_time_c4, v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



input_mat1 = squeeze(struct_in1.gamma_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in4.gamma_data(:,chan3_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 1- Cz gamma';
[gamma_zmap_cz, gamma_zmapthresh_cz, gamma_zmapthresh_plc_cz, gamma_zmapthresh_clc_cz, tfv_time_cz, v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);




% lf
v_freq = struct_in1.lf_frex;
input_mat1 = squeeze(struct_in1.lf_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in4.lf_data(:,chan1_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 1- C3 low frequency';
[lf_zmap_c3, lf_zmapthresh_c3, lf_zmapthresh_plc_c3, lf_zmapthresh_clc_c3, tfv_time_c3, lf_v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

input_mat1 = squeeze(struct_in1.lf_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in4.lf_data(:,chan2_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 1- C4 low frequency';
[lf_zmap_c4, lf_zmapthresh_c4, lf_zmapthresh_plc_c4, lf_zmapthresh_clc_c4, tfv_time_c4, lf_v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


input_mat1 = squeeze(struct_in1.lf_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in4.lf_data(:,chan3_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 1- Cz low frequency';
[lf_zmap_cz, lf_zmapthresh_cz, lf_zmapthresh_plc_cz, lf_zmapthresh_clc_cz, tfv_time_cz, lf_v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



%%%%%%%%%%% Cond 2 vs Cond 3
n_permutes = 1000;
v_time = struct_in1.times;
v_freq = struct_in1.gamma_frex;
time_window = [-150 800];
plotting_input = 1;
input_mat1 = squeeze(struct_in2.gamma_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in3.gamma_data(:,chan1_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 2- C3 gamma';
[gamma_zmap_c3, gamma_zmapthresh_c3, gamma_zmapthresh_plc_c3, gamma_zmapthresh_clc_c3, tfv_time_c3, v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

input_mat1 = squeeze(struct_in2.gamma_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in3.gamma_data(:,chan2_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 2- C4 gamma';
[gamma_zmap_c4, gamma_zmapthresh_c4, gamma_zmapthresh_plc_c4, gamma_zmapthresh_clc_c4, tfv_time_c4, v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



input_mat1 = squeeze(struct_in2.gamma_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in3.gamma_data(:,chan3_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 2- C4 gamma';
[gamma_zmap_cz, gamma_zmapthresh_cz, gamma_zmapthresh_plc_cz, gamma_zmapthresh_clc_cz, tfv_time_cz, v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);




% lf
v_freq = struct_in1.lf_frex;
input_mat1 = squeeze(struct_in2.lf_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in3.lf_data(:,chan1_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 2- C3 low frequency';
[lf_zmap_c3, lf_zmapthresh_c3, lf_zmapthresh_plc_c3, lf_zmapthresh_clc_c3, tfv_time_c3, lf_v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

input_mat1 = squeeze(struct_in2.lf_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in3.lf_data(:,chan2_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 2- C4 low frequency';
[lf_zmap_c4, lf_zmapthresh_c4, lf_zmapthresh_plc_c4, lf_zmapthresh_clc_c4, tfv_time_c4, lf_v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


input_mat1 = squeeze(struct_in2.lf_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in3.lf_data(:,chan3_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 2- C4 low frequency';
[lf_zmap_cz, lf_zmapthresh_cz, lf_zmapthresh_plc_cz, lf_zmapthresh_clc_cz, tfv_time_cz, lf_v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



%%%%%%%%%%% Cond 2 vs Cond 4
n_permutes = 1000;
v_time = struct_in1.times;
v_freq = struct_in1.gamma_frex;
time_window = [-150 800];
plotting_input = 1;
input_mat1 = squeeze(struct_in2.gamma_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in4.gamma_data(:,chan1_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 2- C3 gamma';
[gamma_zmap_c3, gamma_zmapthresh_c3, gamma_zmapthresh_plc_c3, gamma_zmapthresh_clc_c3, tfv_time_c3, v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

input_mat1 = squeeze(struct_in2.gamma_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in4.gamma_data(:,chan2_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 2- C4 gamma';
[gamma_zmap_c4, gamma_zmapthresh_c4, gamma_zmapthresh_plc_c4, gamma_zmapthresh_clc_c4, tfv_time_c4, v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



input_mat1 = squeeze(struct_in2.gamma_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in4.gamma_data(:,chan3_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 2- C4 gamma';
[gamma_zmap_cz, gamma_zmapthresh_cz, gamma_zmapthresh_plc_cz, gamma_zmapthresh_clc_cz, tfv_time_cz, v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);




% lf
v_freq = struct_in1.lf_frex;
input_mat1 = squeeze(struct_in2.lf_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in4.lf_data(:,chan1_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 2- C3 low frequency';
[lf_zmap_c3, lf_zmapthresh_c3, lf_zmapthresh_plc_c3, lf_zmapthresh_clc_c3, tfv_time_c3, lf_v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

input_mat1 = squeeze(struct_in2.lf_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in4.lf_data(:,chan2_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 2- C4 low frequency';
[lf_zmap_c4, lf_zmapthresh_c4, lf_zmapthresh_plc_c4, lf_zmapthresh_clc_c4, tfv_time_c4, lf_v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


input_mat1 = squeeze(struct_in2.lf_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in4.lf_data(:,chan3_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 2- C4 low frequency';
[lf_zmap_cz, lf_zmapthresh_cz, lf_zmapthresh_plc_cz, lf_zmapthresh_clc_cz, tfv_time_cz, lf_v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



%%%%%%%%%%% Cond 3 vs Cond 4
n_permutes = 1000;
v_time = struct_in1.times;
v_freq = struct_in1.gamma_frex;
time_window = [-150 800];
plotting_input = 1;
input_mat1 = squeeze(struct_in3.gamma_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in4.gamma_data(:,chan1_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 3- C3 gamma';
[gamma_zmap_c3, gamma_zmapthresh_c3, gamma_zmapthresh_plc_c3, gamma_zmapthresh_clc_c3, tfv_time_c3, v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

input_mat1 = squeeze(struct_in3.gamma_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in4.gamma_data(:,chan2_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 3- C4 gamma';
[gamma_zmap_c4, gamma_zmapthresh_c4, gamma_zmapthresh_plc_c4, gamma_zmapthresh_clc_c4, tfv_time_c4, v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);



input_mat1 = squeeze(struct_in3.gamma_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in4.gamma_data(:,chan3_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 3- C4 gamma';
[gamma_zmap_cz, gamma_zmapthresh_cz, gamma_zmapthresh_plc_cz, gamma_zmapthresh_clc_cz, tfv_time_cz, v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);




% lf
v_freq = struct_in1.lf_frex;
input_mat1 = squeeze(struct_in3.lf_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in4.lf_data(:,chan1_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 3- C3 low frequency';
[lf_zmap_c3, lf_zmapthresh_c3, lf_zmapthresh_plc_c3, lf_zmapthresh_clc_c3, tfv_time_c3, lf_v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

input_mat1 = squeeze(struct_in3.lf_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in4.lf_data(:,chan2_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 3- C4 low frequency';
[lf_zmap_c4, lf_zmapthresh_c4, lf_zmapthresh_plc_c4, lf_zmapthresh_clc_c4, tfv_time_c4, lf_v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


input_mat1 = squeeze(struct_in3.lf_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in4.lf_data(:,chan3_IDX,:,:));
my_title = 'NPP t-test condition 4 - codition 3- C4 low frequency';
[lf_zmap_cz, lf_zmapthresh_cz, lf_zmapthresh_plc_cz, lf_zmapthresh_clc_cz, tfv_time_cz, lf_v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


%% Doing the plottings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots of signal to noise and difference

% loading data
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data
gamma_in = load('nppc_gamma_new.mat');
lf_in = load('nppc_lf_new.mat');
my_results  = load('results_cond1_cond3.mat');
load('data_grouped_new.mat');
nsubs = 17;

% Gamma - for condition 1 vs 3  - cheps vs leps normal
v_time = cond_1_clean.times;
v_freq = cond_1_clean.gamma_frex;
x = cond_1_clean.times;
time_window = [-150 800];

time_s = dsearchn(v_time',time_window(1)); % start time for tf window
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
x  = x(time_s:time_e);
y = cond_1_clean.gamma_frex;

clims_in = [-5 5];
figure
hold on
subplot(3,3,1)
title_in = ["C3","Cheps"];
zmap = gamma_in.c1_zmap_c3;
mask_in  = gamma_in.c1_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(3,3,2)
title_in = ["Cz","Cheps"];
zmap = gamma_in.c1_zmap_cz;
mask_in  = gamma_in.c1_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(3,3,3)
title_in = ["C4","Cheps"];
zmap = gamma_in.c1_zmap_c4;
mask_in  = gamma_in.c1_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)

subplot(3,3,4)
title_in = ["Leps"];
zmap = gamma_in.c3_zmap_c3;
mask_in  = gamma_in.c3_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(3,3,5)
zmap = gamma_in.c3_zmap_cz;
mask_in  = gamma_in.c3_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(3,3,6)
zmap = gamma_in.c3_zmap_c4;
mask_in  = gamma_in.c3_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
sgtitle('Cheps vs Leps baseline - Gamma')
set(gcf,'color','w');
% Making significant difference plots for gamma for C3, Cz, and C4
my_clims = [-5 5];
my_xlims = [-150 800];
%C3
subplot(3,3,7)
zmap = my_results.gamma_zmap_c3;
mask_in = my_results.gamma_zmapthresh_clc_c3;
my_title = 'Difference (Leps-Cheps)';
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%Cz
subplot(3,3,8)
zmap = my_results.gamma_zmap_cz;
mask_in = my_results.gamma_zmapthresh_clc_cz;
my_title = 'Difference (Leps-Cheps)';
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%C4
subplot(3,3,9)
zmap = my_results.gamma_zmap_c4;
mask_in = my_results.gamma_zmapthresh_clc_c4;
my_title = 'Difference (Leps-Cheps)';
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%colormap(rbw_good_one_2)




%% Low frequency reponse
y = cond_1_clean.lf_frex;
clims_in = [-10 10];

figure
hold on
subplot(3,3,1)
title_in = ["C3","Cheps"];
zmap = lf_in.c1_zmap_c3;
mask_in  = lf_in.c1_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(3,3,2)
title_in = ["Cz","Cheps"];
zmap = lf_in.c1_zmap_cz;
mask_in  = lf_in.c1_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(3,3,3)
title_in = ["C4","Cheps"];
zmap = lf_in.c1_zmap_c4;
mask_in  = lf_in.c1_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)

subplot(3,3,4)
title_in = ["Leps"];
zmap = lf_in.c3_zmap_c3;
mask_in  = lf_in.c3_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(3,3,5)
zmap = lf_in.c3_zmap_cz;
mask_in  = lf_in.c3_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(3,3,6)
zmap = lf_in.c3_zmap_c4;
mask_in  = lf_in.c3_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
sgtitle('Cheps vs Leps baseline - Lf')

% difference
my_clims = [-8 8];
my_xlims = [-150 800];
%C3
subplot(3,3,7)
zmap = my_results.lf_zmap_c3;
mask_in_1 = my_results.lf_zmapthresh_plc_c3;
mask_in_2 = my_results.lf_zmapthresh_clc_c3;
mask_in = mask_in_1+mask_in_2~=0;
my_title = 'Difference (Leps-Cheps)';
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%Cz
subplot(3,3,8)
mask_in_1 = my_results.lf_zmapthresh_plc_cz;
mask_in_2 = my_results.lf_zmapthresh_clc_cz;
mask_in = mask_in_1+mask_in_2~=0;
my_title = 'Difference (Leps-Cheps)';
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%C4
subplot(3,3,9)
mask_in_1 = my_results.lf_zmapthresh_plc_c4;
mask_in_2 = my_results.lf_zmapthresh_clc_c4;
mask_in = mask_in_1+mask_in_2~=0;
mask_in = my_results.lf_zmapthresh_plc_c4;
my_title = 'Difference (Leps-Cheps)';
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%colormap(rbw_good_one_2)
set(gcf,'color','w');


%% overlapping significant signals vs areas of significant differences
%% Gammma
% Cond 1 vs Cond 3 - Cheps vs Leps baseline
cond1 = 'c1';
cond2 = 'c3';
struct_in1 = cond_1_clean;
struct_in2 = cond_3_clean;
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

%%%%%%%%%%%%%%%%%%%%%%%  C3 gamma
mask1 = eval(['gamma_in.',cond1,'_zmapthresh_clc_c3;']);
mask2 = eval(['gamma_in.',cond2,'_zmapthresh_clc_c3;']);
mask_difference = my_results.gamma_zmapthresh_clc_c3~=0;
mask_signal = mask1+mask2~=0;
mask_final = mask_signal+mask_difference==2;

% subset data based on the mask
% getting some parameters
time_window = [-150 800];
v_time = struct_in1.times;

time_s = dsearchn(v_time',time_window(1)); % start time for tf window
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
tfv_time  = v_time(time_s:time_e);
nTimepoints = numel(tfv_time);

if mean(mean(mask_signal)) ~= 0 && mean(mean(mask_difference))~=0
chan_index = 1;
data_in1 = squeeze(struct_in1.gamma_data(:,chan_index,:,time_s:time_e));
data_out1 = zeros(size(data_in1));
data_in2 = squeeze(struct_in2.gamma_data(:,chan_index,:,time_s:time_e));
data_out2 = zeros(size(data_in2));
means_output_cc = zeros(nsubs,2);
figure
for i = 1:17
    s_dat = squeeze(data_in1(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(.001)
    data_out1(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_cc(i,1) = mean(temp);
    
    
    s_dat = squeeze(data_in2(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(.001)
    data_out2(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_cc(i,2) = mean(temp);
end

figure
g = {'Cheps Placebo', 'Leps Placebo'};
boxplot(means_output_cc,g)
my_ylims = [floor(min(min(means_output_cc))-.5) ceil(max(max(means_output_cc))+.5)];
hold on
    for j = 1:nsubs
        y1 = means_output_cc(j,:); % get first condition
        x1 = [1 2];
        lh = plot(x1,y1,'Color',[.05 .05 .05],...
            'LineWidth',0.1,...
            'Marker','.',...
            'MarkerEdgeColor','k',...
            'MarkerSize',12);
        lh.Color =[1,0,0,0.2];
    end
[h,p,ci,stats] = ttest(means_output_cc(:,1),means_output_cc(:,2))
if p<0.05
    yt = get(gca, 'YTick');
    axis([xlim    0  ceil(max(yt)*1.2)])
    xt = get(gca, 'XTick');
    hold on
    plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
end
ylim(my_ylims);
sgtitle('Mean differences in gamma C3')
set(gcf,'color','w');
end





%%%%%%%%%%%%%%%%%% Cz gamma
mask1 = eval(['gamma_in.',cond1,'_zmapthresh_clc_cz;']);
mask2 = eval(['gamma_in.',cond2,'_zmapthresh_clc_cz;']);
mask_difference = my_results.gamma_zmapthresh_clc_cz~=0;
mask_signal = mask1+mask2~=0;
mask_final = mask_signal+mask_difference==2;

if mean(mean(mask_signal)) ~= 0 && mean(mean(mask_difference))~=0
chan_index = 2;
data_in1 = squeeze(struct_in1.gamma_data(:,chan_index,:,time_s:time_e));
data_out1 = zeros(size(data_in1));
data_in2 = squeeze(struct_in2.gamma_data(:,chan_index,:,time_s:time_e));
data_out2 = zeros(size(data_in2));
means_output_cz = zeros(19,2);
figure
for i = 1:nsubs
    s_dat = squeeze(data_in1(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(.1)
    data_out1(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_cz(i,1) = mean(temp);
    
    
    s_dat = squeeze(data_in2(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(.1)
    data_out2(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_cz(i,2) = mean(temp);
end

figure
g = {'Cheps Placebo', 'Leps Placebo'};
my_ylims = [floor(min(min(means_output_cz))-.5) ceil(max(max(means_output_cz))+.5)];
boxplot(means_output_cz,g)
hold on
    for j = 1:nsubs
        y1 = means_output_cz(j,:); % get first condition
        x1 = [1 2];
        lh = plot(x1,y1,'Color',[.05 .05 .05],...
            'LineWidth',0.1,...
            'Marker','.',...
            'MarkerEdgeColor','k',...
            'MarkerSize',12);
        lh.Color =[1,0,0,0.2];
    end
[h,p,ci,stats] = ttest(means_output_cz(:,1),means_output_cz(:,2))
if p<0.05
    yt = get(gca, 'YTick');
    axis([xlim    0  ceil(max(yt)*1.2)])
    xt = get(gca, 'XTick');
    hold on
    plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
end
ylim(my_ylims);
sgtitle('Mean differences in gamma Cz')
set(gcf,'color','w');
end


%%%%%%%%%%%%%%%%%% C4 gamma
mask1 = eval(['gamma_in.',cond1,'_zmapthresh_clc_c4;']);
mask2 = eval(['gamma_in.',cond2,'_zmapthresh_clc_c4;']);
mask_difference = my_results.gamma_zmapthresh_clc_c4~=0;
mask_signal = mask1+mask2~=0;
mask_final = mask_signal+mask_difference==2;

if mean(mean(mask_signal)) ~= 0 && mean(mean(mask_difference))~=0
chan_index = 3;
data_in1 = squeeze(struct_in1.gamma_data(:,chan_index,:,time_s:time_e));
data_out1 = zeros(size(data_in1));
data_in2 = squeeze(struct_in2.gamma_data(:,chan_index,:,time_s:time_e));
data_out2 = zeros(size(data_in2));
means_output_ci = zeros(19,2);
for i = 1:nsubs
    s_dat = squeeze(data_in1(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(.1)
    data_out1(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_ci(i,1) = mean(temp);
    
    
    s_dat = squeeze(data_in2(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(.1)
    data_out2(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_ci(i,2) = mean(temp);
end

figure
g = {'Cheps Placebo', 'Leps Placebo'};
boxplot(means_output_ci,g)
my_ylims = [floor(min(min(means_output_ci))-.5) ceil(max(max(means_output_ci))+.5)];
hold on
    for j = 1:nsubs
        y1 = means_output_ci(j,:); % get first condition
        x1 = [1 2];
        lh = plot(x1,y1,'Color',[.05 .05 .05],...
            'LineWidth',0.1,...
            'Marker','.',...
            'MarkerEdgeColor','k',...
            'MarkerSize',12);
        lh.Color =[1,0,0,0.2];
    end
[h,p,ci,stats] = ttest(means_output_ci(:,1),means_output_ci(:,2))
if p<0.05
    yt = get(gca, 'YTick');
    axis([xlim    0  ceil(max(yt)*1.2)])
    xt = get(gca, 'XTick');
    hold on
    plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
end
ylim(my_ylims);
sgtitle('Mean differences in gamma Cz')
set(gcf,'color','w');
end


%%%%%%%%%%%%%% Plotting comparisons LF %%%%%%%%%%%%%%%
% Lf
% Cond 1 vs Cond 3 - Cheps vs Leps baseline

%%%%%%%%%%%%%%%%%%%%%%%  C3 lf
mask1 = eval(['lf_in.',cond1,'_zmapthresh_clc_c3;']);
mask2 = eval(['lf_in.',cond2,'_zmapthresh_clc_c3;']);
mask_dif_1 = my_results.lf_zmapthresh_plc_c3;
mask_dif_2 = my_results.lf_zmapthresh_clc_c3;
mask_difference = mask_dif_1+mask_dif_2~=0;
%mask_difference = my_results.lf_zmapthresh_clc_c3~=0;
mask_signal = mask1+mask2~=0;
mask_final = mask_signal+mask_difference==2;

% subset data based on the mask
% getting some parameters
time_window = [-150 800];
v_time = struct_in1.times;

time_s = dsearchn(v_time',time_window(1)); % start time for tf window
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
tfv_time  = v_time(time_s:time_e);
nTimepoints = numel(tfv_time);

if mean(mean(mask_signal)) ~= 0 && mean(mean(mask_difference))~=0
chan_index = 1;
data_in1 = squeeze(struct_in1.lf_data(:,chan_index,:,time_s:time_e));
data_out1 = zeros(size(data_in1));
data_in2 = squeeze(struct_in2.lf_data(:,chan_index,:,time_s:time_e));
data_out2 = zeros(size(data_in2));
means_output_cc = zeros(nsubs,2);
for i = 1:nsubs
    s_dat = squeeze(data_in1(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(.001)
    data_out1(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_cc(i,1) = mean(temp);
    
    
    s_dat = squeeze(data_in2(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(.001)
    data_out2(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_cc(i,2) = mean(temp);
end

figure
g = {'Cheps Placebo', 'Leps Placebo'};
boxplot(means_output_cc,g)
my_ylims = [floor(min(min(means_output_cc))-.5) ceil(max(max(means_output_cc))+.5)];
hold on
    for j = 1:nsubs
        y1 = means_output_cc(j,:); % get first condition
        x1 = [1 2];
        lh = plot(x1,y1,'Color',[.05 .05 .05],...
            'LineWidth',0.1,...
            'Marker','.',...
            'MarkerEdgeColor','k',...
            'MarkerSize',12);
        lh.Color =[1,0,0,0.2];
    end
[h,p,ci,stats] = ttest(means_output_cc(:,1),means_output_cc(:,2))
if p<0.05
    yt = get(gca, 'YTick');
    axis([xlim    0  ceil(max(yt)*1.2)])
    xt = get(gca, 'XTick');
    hold on
    plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
end
ylim(my_ylims);
sgtitle('Mean differences in Lf C3')
end





%%%%%%%%%%%%%%%%%% Cz Lf
mask1 = eval(['lf_in.',cond1,'_zmapthresh_clc_cz;']);
mask2 = eval(['lf_in.',cond2,'_zmapthresh_clc_cz;']);
mask_dif_1 = my_results.lf_zmapthresh_plc_cz;
mask_dif_2 = my_results.lf_zmapthresh_clc_cz;
mask_difference = mask_dif_1+mask_dif_2~=0;
%mask_difference = my_results.lf_zmapthresh_clc_cz~=0;
mask_signal = mask1+mask2~=0;
mask_final = mask_signal+mask_difference==2;

if mean(mean(mask_signal)) ~= 0 && mean(mean(mask_difference))~=0
chan_index = 2;
data_in1 = squeeze(struct_in1.lf_data(:,chan_index,:,time_s:time_e));
data_out1 = zeros(size(data_in1));
data_in2 = squeeze(struct_in2.lf_data(:,chan_index,:,time_s:time_e));

data_out2 = zeros(size(data_in2));
means_output_cz = zeros(19,2);
figure
for i = 1:nsubs
    s_dat = squeeze(data_in1(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(.1)
    data_out1(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_cz(i,1) = mean(temp);
    
    
    s_dat = squeeze(data_in2(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(.1)
    data_out2(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_cz(i,2) = mean(temp);
end

figure
g = {'Cheps Placebo', 'Leps Placebo'};
my_ylims = [floor(min(min(means_output_cz))-.5) ceil(max(max(means_output_cz))+.5)];
boxplot(means_output_cz,g)
hold on
    for j = 1:nsubs
        y1 = means_output_cz(j,:); % get first condition
        x1 = [1 2];
        lh = plot(x1,y1,'Color',[.05 .05 .05],...
            'LineWidth',0.1,...
            'Marker','.',...
            'MarkerEdgeColor','k',...
            'MarkerSize',12);
        lh.Color =[1,0,0,0.2];
    end
[h,p,ci,stats] = ttest(means_output_cz(:,1),means_output_cz(:,2))
if p<0.05
    yt = get(gca, 'YTick');
    axis([xlim    0  ceil(max(yt)*1.2)])
    xt = get(gca, 'XTick');
    hold on
    plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
end
ylim(my_ylims);
sgtitle('Mean differences in Lf Cz')
end


%%%%%%%%%%%%%%%%%% C4 lf
mask1 = eval(['lf_in.',cond1,'_zmapthresh_clc_c4;']);
mask2 = eval(['lf_in.',cond2,'_zmapthresh_clc_c4;']);
%mask_difference = my_results.gamma_zmapthresh_clc_c4~=0;
mask_dif_1 = my_results.lf_zmapthresh_plc_c4;
mask_dif_2 = my_results.lf_zmapthresh_clc_c4;
mask_difference = mask_dif_1+mask_dif_2~=0;
mask_signal = mask1+mask2~=0;
mask_final = mask_signal+mask_difference==2;

if mean(mean(mask_signal)) ~= 0 && mean(mean(mask_difference))~=0
chan_index = 3;
data_in1 = squeeze(struct_in1.lf_data(:,chan_index,:,time_s:time_e));
data_out1 = zeros(size(data_in1));
data_in2 = squeeze(struct_in2.lf_data(:,chan_index,:,time_s:time_e));
data_out2 = zeros(size(data_in2));
means_output_ci = zeros(19,2);
for i = 1:nsubs
    s_dat = squeeze(data_in1(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(.1)
    data_out1(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_ci(i,1) = mean(temp);
    
    
    s_dat = squeeze(data_in2(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(.1)
    data_out2(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_ci(i,2) = mean(temp);
end

figure
g = {'Cheps Placebo', 'Leps Placebo'};
boxplot(means_output_ci,g)
my_ylims = [floor(min(min(means_output_ci))-.5) ceil(max(max(means_output_ci))+.5)];
hold on
    for j = 1:nsubs
        y1 = means_output_ci(j,:); % get first condition
        x1 = [1 2];
        lh = plot(x1,y1,'Color',[.05 .05 .05],...
            'LineWidth',0.1,...
            'Marker','.',...
            'MarkerEdgeColor','k',...
            'MarkerSize',12);
        lh.Color =[1,0,0,0.2];
    end
[h,p,ci,stats] = ttest(means_output_ci(:,1),means_output_ci(:,2))
if p<0.05
    yt = get(gca, 'YTick');
    axis([xlim    0  ceil(max(yt)*1.2)])
    xt = get(gca, 'XTick');
    hold on
    plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
end
ylim(my_ylims);
sgtitle('Mean differences in lf C4')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cheps Placebo vs cheps caps
my_results  = load('results_cond1_cond2.mat');


% Gamma - for condition 1 vs 2  - cheps vs  cheps caps
x = cond_1_clean.times;
time_s = dsearchn(v_time',time_window(1)); % start time for tf window
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
x  = x(time_s:time_e);
y = cond_1_clean.gamma_frex;

clims_in = [-5 5];
figure
hold on
subplot(3,3,1)
title_in = ["C3","Cheps Baseline"];
zmap = gamma_in.c1_zmap_c3;
mask_in  = gamma_in.c1_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(3,3,2)
title_in = ["Cz","Cheps Baseline"];
zmap = gamma_in.c1_zmap_cz;
mask_in  = gamma_in.c1_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(3,3,3)
title_in = ["C4","Cheps Baseline"];
zmap = gamma_in.c1_zmap_c4;
mask_in  = gamma_in.c1_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)

subplot(3,3,4)
title_in = ["Cheps Cap"];
zmap = gamma_in.c2_zmap_c3;
mask_in  = gamma_in.c2_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(3,3,5)
zmap = gamma_in.c2_zmap_cz;
mask_in  = gamma_in.c2_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(3,3,6)
zmap = gamma_in.c2_zmap_c4;
mask_in  = gamma_in.c2_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
sgtitle('Cheps Baseline vs Cheps cap - Gamma')

% Making significant difference plots for gamma for C3, Cz, and C4
my_clims = [-5 5];
my_xlims = [-150 800];
%C3
subplot(3,3,7)
zmap = my_results.gamma_zmap_c3;
mask_in = my_results.gamma_zmapthresh_clc_c3;
my_title = 'Difference (Cheps caps - Cheps baseline)';
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%Cz
subplot(3,3,8)
zmap = my_results.gamma_zmap_cz;
mask_in = my_results.gamma_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%C4
subplot(3,3,9)
zmap = my_results.gamma_zmap_c4;
mask_in = my_results.gamma_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%colormap(rbw_good_one_2)




% Low frequency reponse
y = cond_1_clean.lf_frex;
clims_in = [-10 10];

figure
hold on
subplot(3,3,1)
title_in = ["C3","Cheps Baseline"];
zmap = lf_in.c1_zmap_c3;
mask_in  = lf_in.c1_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(3,3,2)
title_in = ["Cz","Cheps Baseline"];
zmap = lf_in.c1_zmap_cz;
mask_in  = lf_in.c1_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(3,3,3)
title_in = ["C4","Cheps Baseline"];
zmap = lf_in.c1_zmap_c4;
mask_in  = lf_in.c1_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)

subplot(3,3,4)
title_in = ["Cheps Cap"];
zmap = lf_in.c2_zmap_c3;
mask_in  = lf_in.c2_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(3,3,5)
zmap = lf_in.c2_zmap_cz;
mask_in  = lf_in.c2_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(3,3,6)
zmap = lf_in.c2_zmap_c4;
mask_in  = lf_in.c2_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
sgtitle('Cheps Baseline vs Cheps cap - Lf')

% difference
my_clims = [-8 8];
my_xlims = [-150 800];
%C3
subplot(3,3,7)
zmap = my_results.lf_zmap_c3;
mask_in_1 = my_results.lf_zmapthresh_plc_c3;
mask_in_2 = my_results.lf_zmapthresh_clc_c3;
mask_in = mask_in_1+mask_in_2~=0;
my_title = 'Difference (Cheps caps - Cheps baseline)';
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%Cz
subplot(3,3,8)
mask_in_1 = my_results.lf_zmapthresh_plc_cz;
mask_in_2 = my_results.lf_zmapthresh_clc_cz;
mask_in = mask_in_1+mask_in_2~=0;
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%C4
subplot(3,3,9)
mask_in_1 = my_results.lf_zmapthresh_plc_c4;
mask_in_2 = my_results.lf_zmapthresh_clc_c4;
mask_in = mask_in_1+mask_in_2~=0;
mask_in = my_results.lf_zmapthresh_plc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%colormap(rbw_good_one_2)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Leps Placebo vs Leps caps
my_results  = load('results_cond3_cond4.mat');


% Gamma - for condition 1 vs 2  - cheps vs  cheps caps
x = cond_1_clean.times;
time_s = dsearchn(v_time',time_window(1)); % start time for tf window
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
x  = x(time_s:time_e);
y = cond_1_clean.gamma_frex;

clims_in = [-5 5];
figure
hold on
subplot(3,3,1)
title_in = ["C3","Leps Baseline"];
zmap = gamma_in.c3_zmap_c3;
mask_in  = gamma_in.c3_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(3,3,2)
title_in = ["Cz","Leps Baseline"];
zmap = gamma_in.c3_zmap_cz;
mask_in  = gamma_in.c3_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(3,3,3)
title_in = ["C4","Leps Baseline"];
zmap = gamma_in.c3_zmap_c4;
mask_in  = gamma_in.c3_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)

subplot(3,3,4)
title_in = ["Leps Cap"];
zmap = gamma_in.c4_zmap_c3;
mask_in  = gamma_in.c4_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(3,3,5)
zmap = gamma_in.c4_zmap_cz;
mask_in  = gamma_in.c4_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(3,3,6)
zmap = gamma_in.c4_zmap_c4;
mask_in  = gamma_in.c4_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
sgtitle('Leps Baseline vs Leps cap - Gamma')

% Making significant difference plots for gamma for C3, Cz, and C4
my_clims = [-5 5];
my_xlims = [-150 800];
%C3
subplot(3,3,7)
zmap = my_results.gamma_zmap_c3;
mask_in = my_results.gamma_zmapthresh_clc_c3;
my_title = 'Difference (Leps caps - Leps baseline)';
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%Cz
subplot(3,3,8)
zmap = my_results.gamma_zmap_cz;
mask_in = my_results.gamma_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%C4
subplot(3,3,9)
zmap = my_results.gamma_zmap_c4;
mask_in = my_results.gamma_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%colormap(rbw_good_one_2)




% Low frequency reponse
y = cond_1_clean.lf_frex;
clims_in = [-10 10];

figure
hold on
subplot(3,3,1)
title_in = ["C3","Leps Baseline"];
zmap = lf_in.c3_zmap_c3;
mask_in  = lf_in.c3_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(3,3,2)
title_in = ["Cz","Leps Baseline"];
zmap = lf_in.c3_zmap_cz;
mask_in  = lf_in.c3_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(3,3,3)
title_in = ["C4","Leps Baseline"];
zmap = lf_in.c3_zmap_c4;
mask_in  = lf_in.c3_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)

subplot(3,3,4)
title_in = ["Leps Cap"];
zmap = lf_in.c4_zmap_c3;
mask_in  = lf_in.c4_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(3,3,5)
zmap = lf_in.c4_zmap_cz;
mask_in  = lf_in.c4_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(3,3,6)
zmap = lf_in.c4_zmap_c4;
mask_in  = lf_in.c4_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
sgtitle('Leps Baseline vs Leps cap - Lf')

% difference
my_clims = [-8 8];
my_xlims = [-150 800];
%C3
subplot(3,3,7)
zmap = my_results.lf_zmap_c3;
mask_in_1 = my_results.lf_zmapthresh_plc_c3;
mask_in_2 = my_results.lf_zmapthresh_clc_c3;
mask_in = mask_in_1+mask_in_2~=0;
my_title = 'Difference (Leps caps - Leps baseline)';
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%Cz
subplot(3,3,8)
mask_in_1 = my_results.lf_zmapthresh_plc_cz;
mask_in_2 = my_results.lf_zmapthresh_clc_cz;
mask_in = mask_in_1+mask_in_2~=0;
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%C4
subplot(3,3,9)
mask_in_1 = my_results.lf_zmapthresh_plc_c4;
mask_in_2 = my_results.lf_zmapthresh_clc_c4;
mask_in = mask_in_1+mask_in_2~=0;
mask_in = my_results.lf_zmapthresh_plc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,my_clims,my_title)
%colormap(rbw_good_one_2)

%% 
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\github_repo
[map,num,typ,scheme] = brewermap_view

my_lims = [-40 40];
x = repmat(EEG_Mast_all_finalized.times',[1, 19]);
y = repmat([1:19],1500,1);
z =  squeeze(EEG_Mast_all_finalized.data_subs(1:19,23,:))';
figure
hold on
for i  = 1:19
   plot3(y(:,i),x(:,i),z(:,i),'Color',map(i,:),'LineWidth',2);
end
ylabel('time (ms)');
xlabel('subject')
zlabel('amplitude (micro volts');
zlim(my_lims);
set(gcf,'color','w');
view(45,-45)
exportgraphics(gcf,'leps_time.eps','ContentType','vector')


x = repmat(EEG_Mast_all_finalized.times',[1, 19]);
y = repmat([1:19],1500,1);
z =  squeeze(EEG_Mast_all_finalized.data_subs(1:19,24,:))';
figure
hold on
for i  = 1:19
   plot3(y(:,i),x(:,i),z(:,i),'Color',map(i,:),'LineWidth',2);
end
ylabel('time (ms)');
xlabel('subject')
zlabel('amplitude (micro volts');
zlim(my_lims);
set(gcf,'color','w');
view(45,-45)
exportgraphics(gcf,'cheps_time.eps','ContentType','vector')



