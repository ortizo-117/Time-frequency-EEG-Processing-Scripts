%% statistics bruv

cd Z:\30_Oscar_Ortiz\GluEEG\gamma_sta
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\scripts_functions



%% Determing ROI of activity

interp_fact_x = 4;
interp_fact_y = 2;
time_window = [-400 900];
n_permutes = 10;
plotting_window = 1;


struct_in = EEG_Fz_all;
struct_in = interp_matrix_eeg(struct_in,interp_fact_x,interp_fact_y);

EEG_Fz_all_results = wrap_npp_stb_eegstruct(struct_in,n_permutes,time_window,plotting_window);

save('results_Fz','EEG_Fz_all_results')
clear EEG_Fz_all_results EEG_Fz_all


struct_in = EEG_Mast_all;
struct_in = interp_matrix_eeg(struct_in,interp_fact_x,interp_fact_y);
plotting_input = 1;
EEG_Mast_all_results = wrap_npp_stb_eegstruct(struct_in,n_permutes,time_window,plotting_input);

save('results_Mast','EEG_Mast_all_results')
clear EEG_Mast_all_results EEG_Mast_all


struct_in = EEG_AVG_all;
struct_in = interp_matrix_eeg(struct_in,interp_fact_x,interp_fact_y);
plotting_input = 1;
EEG_AVG_all_results = wrap_npp_stb_eegstruct(struct_in,n_permutes,time_window,plotting_input);

save('results_AVG.mat','EEG_AVG_all_results')
clear EEG_AVG_all_results EEG_AVG_all


%% Plotting gamma and significant ROIs
my_chanels = EEG_Mast_all_finalized.chanlocs;


chan_name1 = "C3";
chan_name2 = "Cz";
chan_name3 = "C4";
chan_name4 = "Fc1";
chan_name5 = "Fc2";
chan_name6 = "Pz";
chan_name7 = "O2";
chan_name8 = "O1";

for k = 1:length(my_chanels)
    C_chan = lower(string(my_chanels(k).labels));
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
    if C_chan == lower(chan_name6)
        chan6_IDX=k;
    end
    if C_chan == lower(chan_name7)
        chan7_IDX=k;
    end
    if C_chan == lower(chan_name8)
        chan8_IDX=k;
    end
end

struct_in = EEG_Mast_all_finalized;
clear EEG_Mast_all_results


my_clims = [-1.6 1.6];
my_xlims = [-200 800];
x = struct_in.tf_time;
x1 = struct_in.zmaps.time_v;
y = struct_in.tf_frex;


figure 
subplot(4,3,4);
data_in = squeeze(struct_in.tf_data(chan1_IDX,:,:));
mask_in = squeeze(struct_in.zmaps.zmapthresh_clc(chan1_IDX,:,:));
title_in = struct_in.chanlocs(chan1_IDX).labels;
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,my_clims,my_xlims,title_in)

subplot(4,3,5);
data_in = squeeze(struct_in.tf_data(chan2_IDX,:,:));
mask_in = squeeze(struct_in.zmaps.zmapthresh_clc(chan2_IDX,:,:));
title_in = struct_in.chanlocs(chan2_IDX).labels;
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,my_clims,my_xlims,title_in)


subplot(4,3,6);
data_in = squeeze(struct_in.tf_data(chan3_IDX,:,:));
mask_in = squeeze(struct_in.zmaps.zmapthresh_clc(chan3_IDX,:,:));
title_in = struct_in.chanlocs(chan3_IDX).labels;
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,my_clims,my_xlims,title_in)


% subplot(3,3,2);
% data_in = squeeze(struct_in.tf_data(chan4_IDX,:,:));
% mask_in = squeeze(struct_in.zmaps.zmapthresh_clc(chan4_IDX,:,:));
% title_in = struct_in.chanlocs(chan4_IDX).labels;
% plot_og_dat_with_mask(data_in,mask_in,x,x1,y,my_clims,my_xlims,title_in)


subplot(4,3,8);
data_in = squeeze(struct_in.tf_data(chan6_IDX,:,:));
mask_in = squeeze(struct_in.zmaps.zmapthresh_clc(chan6_IDX,:,:));
title_in = struct_in.chanlocs(chan6_IDX).labels;
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,my_clims,my_xlims,title_in)


subplot(4,3,1);
data_in = squeeze(struct_in.tf_data(chan4_IDX,:,:));
mask_in = squeeze(struct_in.zmaps.zmapthresh_clc(chan4_IDX,:,:));
title_in = struct_in.chanlocs(chan4_IDX).labels;
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,my_clims,my_xlims,title_in)


subplot(4,3,3);
data_in = squeeze(struct_in.tf_data(chan5_IDX,:,:));
mask_in = squeeze(struct_in.zmaps.zmapthresh_clc(chan5_IDX,:,:));
title_in = struct_in.chanlocs(chan5_IDX).labels;
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,my_clims,my_xlims,title_in)

subplot(4,3,10);
data_in = squeeze(struct_in.tf_data(chan8_IDX,:,:));
mask_in = squeeze(struct_in.zmaps.zmapthresh_clc(chan8_IDX,:,:));
title_in = struct_in.chanlocs(chan8_IDX).labels;
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,my_clims,my_xlims,title_in)


subplot(4,3,12);
data_in = squeeze(struct_in.tf_data(chan7_IDX,:,:));
mask_in = squeeze(struct_in.zmaps.zmapthresh_clc(chan7_IDX,:,:));
title_in = struct_in.chanlocs(chan7_IDX).labels;
plot_og_dat_with_mask(data_in,mask_in,x,x1,y,my_clims,my_xlims,title_in)

set(gcf,'color','white')





%% Extracting 


cd Z:\30_Oscar_Ortiz\GluEEG\gamma_sta
load('results_Mast.mat')
time_windows_to_extract = [200];
length_window = 200;
EEG_Mast_all_finalized = erd_ers_extraction_struct(EEG_Mast_all_results,time_windows_to_extract,length_window);


%% Topoplotting
addpath C:\Users\ortizo\Documents\MATLAB\eeglab2021.1
eeglab
cd Z:\30_Oscar_Ortiz\GluEEG\gamma_sta
%load("results_Mast_final.mat");

% weighted ERS
map_lims =[-.35 .35];
values_in = squeeze(mean(EEG_Mast_all_finalized.extraction_results.w_mean_ERS));
values_in(1,:) = [];
loc_file = EEG_Mast_all_finalized.chanlocs;
loc_file(1) = [];
n_windows = length(EEG_Mast_all_finalized.extraction_results.time_windows);
figure
sgtitle("weighted ers")
for i = 1:n_windows
    ers_toplot = values_in(:,i);
    subplot(1,n_windows,i)
    topoplot(ers_toplot,loc_file,'maplimits',map_lims);
    title(num2str(EEG_Mast_all_finalized.extraction_results.time_windows(i)))
end



% weighted ERD
map_lims =[-.1 .1];
values_in = squeeze(mean(EEG_Mast_all_finalized.extraction_results.w_mean_ERD));
values_in(1,:) = [];
loc_file = EEG_Mast_all_finalized.chanlocs;
loc_file(1) = [];
n_windows = length(EEG_Mast_all_finalized.extraction_results.time_windows);
figure
sgtitle("weighted erd")
for i = 1:n_windows
    ers_toplot = values_in(:,i);
    subplot(1,n_windows,i)
    topoplot(ers_toplot,loc_file,'maplimits',map_lims);
    title(num2str(EEG_Mast_all_finalized.extraction_results.time_windows(i)))
end



% Mean ERS
map_lims =[-.2 .6];
values_in = squeeze(mean(EEG_Mast_all_finalized.extraction_results.mean_ERS));
values_in(1,:) = [];
loc_file = EEG_Mast_all_finalized.chanlocs;
loc_file(1) = [];
n_windows = length(EEG_Mast_all_finalized.extraction_results.time_windows);
figure
sgtitle("mean ERS")
for i = 1:n_windows
    ers_toplot = values_in(:,i);
    subplot(1,n_windows,i)
    %topoplot(ers_toplot,loc_file,'maplimits',map_lims);
    headplot(ers_toplot, 'STUDY_headplot.spl', headplotparams{:}, 'maplimits', map_lims, 'lighting', 'on');
    title(num2str(EEG_Mast_all_finalized.extraction_results.time_windows(i)))
end

% Mean ERD
map_lims =[-.4 .1];
values_in = squeeze(mean(EEG_Mast_all_finalized.extraction_results.mean_ERD));
values_in(1,:) = [];
loc_file = EEG_Mast_all_finalized.chanlocs;
loc_file(1) = [];
n_windows = length(EEG_Mast_all_finalized.extraction_results.time_windows);
figure
sgtitle("mean ERD")
for i = 1:n_windows
    ers_toplot = values_in(:,i);
    subplot(1,n_windows,i)
    topoplot(ers_toplot,loc_file,'maplimits',map_lims);
    title(num2str(EEG_Mast_all_finalized.extraction_results.time_windows(i)))
end


%% 3d headplots

%% Simple 3-D movie
% Use the graphic interface to coregister your head model with your electrode positions
headplotparams1 = { 'meshfile', 'mheadnew.mat'       , 'transform', [0.664455     -3.39403     -14.2521  -0.00241453     0.015519     -1.55584           11      10.1455           12] };
headplotparams2 = { 'meshfile', 'colin27headmesh.mat', 'transform', [0          -13            0          0.1            0        -1.57         11.7         12.5           12] };
headplotparams  = headplotparams1; % switch here between 1 and 2

% set up the spline file
headplot('setup', loc_file, 'STUDY_headplot.spl', headplotparams{:}); close
 
figure

my_subs = {'P01','S01','S02','S04','S05','S06','S07','S08','S09','S11','S12','S15','S16','S17','S18','S19','S20','S21','S22','S23','S24','S25','S26','S28','S29'};


%% plotting n2 p2s
my_hex_color = '#05AEA2';
my_color = sscanf(my_hex_color(2:end),'%2x%2x%2x',[1 3])/255;
my_lims_N2P2 = [-25 25];
figure
my_times_ERP = EEG_Mast_all_finalized.times;
N2_P2_data_pre =squeeze(EEG_Mast_all_finalized.data_subs(:,chan2_IDX,:));
N2P2_pre = mean(N2_P2_data_pre);
SEM=std(N2_P2_data_pre)/sqrt(size(N2_P2_data_pre,1));    %Standard Error
CI95=tinv([0.025 0.975],size(N2_P2_data_pre,1)-1);                               %T-Score
ERP_CI95=bsxfun(@times, SEM,CI95(:));  %Confidence Intervals 
lower_pre=N2P2_pre+ERP_CI95(1,:);
upper_pre=N2P2_pre+ERP_CI95(2,:);
plot(my_times_ERP,N2P2_pre,'color','k','LineWidth', 1.5); hold on;
plot(my_times_ERP,upper_pre,'color','w');
hold on; 
plot(my_times_ERP,lower_pre,'color','w'); 
x2=[my_times_ERP,fliplr(my_times_ERP)];
inbetween=[upper_pre, fliplr(lower_pre)];
fill(x2,inbetween,my_color,'FaceAlpha',0.4,'EdgeColor','none'); %set colour black, transparency 10%, and no edgecolour
plot(my_times_ERP,N2P2_pre,'color','k','LineWidth', 1.5); hold on;
ylim(my_lims_N2P2);
%xlim(my_xlims);
set(gca, 'YDir','reverse')
xlabel('Times (ms)')
ylabel('\muV')
set(gcf,'color','w');

line([0 0], [-100 101], 'Color','black','LineStyle','--')


