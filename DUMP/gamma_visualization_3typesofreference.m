%% Load the data per subjet
load("S08_tf_data__FZ_ref.mat")
EEG_Fz = EEG;
load("S08_tf_data__Avg.mat")
EEG_AVG = EEG;
load("S08_tf_data__Mast.mat")
EEG_MAST = EEG;
S_ID = '08';




%% Loading all subject data
%% no data = S10, S13, 14

my_subs = {'P01','S01','S02','S04','S05','S06','S07','S08','S09','S11','S12','S15','S16','S17','S18','S19','S20','S21','S22'};
chans = ["C3", "Cz", "C4"];
my_xlims = [-400 920];
my_clims = [-4 4];
n_freq = 43;
n_time = 1000;
n_time_time = 1500;
n_chan = 31;
EEG_Fz_all_d = zeros([length(my_subs), n_chan, n_freq,n_time]);
EEG_AVG_all_d = EEG_Fz_all_d;
EEG_Mast_all_d = EEG_Fz_all_d;
EEG_Fz_all_time_d = zeros([length(my_subs), n_chan,n_time_time]);
EEG_AVG_all_time_d =EEG_Fz_all_time_d;
EEG_Mast_all_time_d =EEG_Fz_all_time_d;



for i = 1:length(my_subs)
    S_ID = char(my_subs(i));
    load([S_ID '_tf_data__FZ_ref.mat'])
    EEG_Fz = EEG;
    load([S_ID '_tf_data__Avg.mat'])
    EEG_AVG = EEG;
    load([S_ID '_tf_data__Mast.mat'])
    EEG_MAST = EEG;

    EEG_Fz_all_d(i,:,:,:) = EEG_Fz.tf_data;
    EEG_AVG_all_d(i,:,:,:) = EEG_AVG.tf_data;
    EEG_Mast_all_d(i,:,:,:) = EEG_MAST.tf_data;

    EEG_Fz_all_time_d(i,:,:) = squeeze(mean(EEG_Fz.data,3));
    EEG_AVG_all_time_d(i,:,:) = squeeze(mean(EEG_AVG.data,3));
    EEG_Mast_all_time_d(i,:,:) = squeeze(mean(EEG_MAST.data,3));

    
    plot_single_subject_diff_ref(EEG_Fz,EEG_AVG,EEG_MAST,chans,S_ID,my_xlims,my_clims)
end

EEG_Fz_all_mean = squeeze(mean(EEG_Fz_all_d));
EEG_AVG_all_mean = squeeze(mean(EEG_AVG_all_d));
EEG_Mast_all_mean = squeeze(mean(EEG_Mast_all_d));


EEG_Fz_all = EEG_Fz;
EEG_Fz_all.tf_data = EEG_Fz_all_mean;
EEG_Fz_all.tf_data_subs = EEG_Fz_all_d;
EEG_Fz_all.data = squeeze(mean(EEG_Fz_all_time_d));
EEG_Fz_all.data_subs = EEG_Fz_all_time_d;

EEG_AVG_all = EEG_AVG;
EEG_AVG_all.tf_data = EEG_AVG_all_mean;
EEG_AVG_all.tf_data_subs = EEG_AVG_all_d;
EEG_AVG_all.data = squeeze(mean(EEG_AVG_all_time_d));
EEG_AVG_all.data_subs = EEG_AVG_all_time_d;



EEG_Mast_all = EEG_MAST;
EEG_Mast_all.tf_data = EEG_Mast_all_mean;
EEG_Mast_all.tf_data_subs = EEG_Mast_all_d;
EEG_Mast_all.data = squeeze(mean(EEG_Mast_all_time_d));
EEG_Mast_all.data_subs = EEG_Mast_all_time_d;


clear EEG_Mast_all_d EEG_AVG_all_d EEG_Fz_all_d




%% pLOTTING

figure
subplot(3,1,1)
plot(EEG_Fz_all.times,EEG_Fz_all.data')
title('Fz referenced')
subplot(3,1,2)
plot(EEG_Fz_all.times,EEG_AVG_all.data')
title('AVG referenced')
subplot(3,1,3)
plot(EEG_Fz_all.times,EEG_Mast_all.data')
title('Mast referenced')

% Plotting time frequency
my_subs = {'P01','S01','S02','S04','S05','S06','S07','S08','S09','S11','S12','S15','S16','S17','S18','S19','S20','S21','S22'};
chans = ["C3", "Cz", "C4"];
my_xlims = [-400 920];
my_clims = [-2 2];

    

plot_single_subject_diff_ref(EEG_Fz_all,EEG_AVG_all,EEG_Mast_all,chans,S_ID,my_xlims,my_clims)



