%% Author: Oscar Ortiz & Cassie Choles
%  Date: 10/07/2021
%  Purpose: Pre-process LEPS vs. CHEPS gamma pipeline


%% Indicate path with all folders 
clear all 
close all

% addpath C:\Users\Kip\Documents\eeglab2021.1 % add eeglab path
cd  Z:\30_Oscar_Ortiz\LEPs_CHEPs\scripts_functions
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\Raw
addpath C:\Users\User\OneDrive\Documents\MATLAB\eeglab2021.0 % for oscar at home
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\Raw\Participant_C06_Trial 1,3 & cond 2, 4'; % folder 
addpath C:\Users\Kip\Documents\eeglab2021.1
addpath C:\Users\ortizo\Documents\MATLAB\eeglab2021.1
%% preprocess
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\0_raw_data\actual_data';
saving_folder = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\1_preprocessed_data\round2';
extension_to_look_for = '*.set';
preproces_EEG_folder(pathname, saving_folder,extension_to_look_for)

%% blink and artifact detection step
addpath C:\Users\Kip\Documents\eeglab2021.1;
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\1_preprocessed_data\round2';
saving_folder = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\2_cleaned_data\round2';
%manual_blink_rejection(pathname);

file_struct_list = dir([pathname filesep() '*ref_filt_ICA_epoch.set']);  %% get list of .set files in the pathname specified

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

pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\2_cleaned_data\round2\already_processed_original'; % folder 

%files_not_processed = gamma_detection(pathname);
epoch_length = [-500 1000];
freqs_in = [1 100];
baseline_in = [-450 -50];
n_times_out = 1000;
pads = 4 ;% can be 2^n 
win_size = 150;
NF = [];
c_lims = [-4 4];

fiels_not_processed = gamma_detection_v2(pathname,epoch_length,freqs_in,baseline_in,n_times_out,pads,win_size,NF,c_lims)

%% Grouping of time-frequency data
addpath Z:\29_Hannah_Goodings\1_Projects\CAP_LONG\CHEP_Only\3_Data_Export
my_chans_in = {'c3','cz','c4','fz','pz'};
% Condition 1
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\round2\Cond1';
cond_1 = grouping_data_specify_chans(pathname,my_chans_in);

% Condition 2
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\round2\Cond2';
cond_2 =  grouping_data_specify_chans(pathname,my_chans_in);

% Condition 3
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\round2\Cond3';
cond_3 =  grouping_data_specify_chans(pathname,my_chans_in);

% Condition 4
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\round2\Cond4';
cond_4 =  grouping_data_specify_chans(pathname,my_chans_in);

%% Grouping N2_P2 data
addpath Z:\29_Hannah_Goodings\1_Projects\CAP_LONG\CHEP_Only\3_Data_Export
addpath C:\Users\ortizo\Documents\MATLAB\eeglab2021.1

my_chans_in = {'c3','cz','c4','fz','pz'};
% Condition 1
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\round2\Cond1';
cond_1 = grouping_data_specify_chans_evoked_pot(pathname,my_chans_in);

% Condition 2
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\round2\Cond2';
cond_2 =  grouping_data_specify_chans_evoked_pot(pathname,my_chans_in);

% Condition 3
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\round2\Cond3';
cond_3 =  grouping_data_specify_chans_evoked_pot(pathname,my_chans_in);

% Condition 4
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\round2\Cond4';
cond_4 =  grouping_data_specify_chans_evoked_pot(pathname,my_chans_in);



%% Plotting_n2_p2
% Creating pulses to plot over  
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data\round2

load('n2_p2_data.mat')

my_time = cond_1.times;
% cheps

% lowpass filter the data

figure
my_chans_in = {'c3','cz','c4','fz','pz'};
my_lims_N2P2 = [-300 800];
my_xlims = [-22 22];
chan_idx = 2;
N2_P2 =squeeze(cond_1.data(:,chan_idx,:));
mean_N2_P2 = mean(N2_P2);
plot(my_time,mean_N2_P2,'color','b','LineWidth', 1.5); hold on;
N2_P2 =squeeze(cond_3.data(:,chan_idx,:));
mean_N2_P2 = mean(N2_P2);
plot(my_time,mean_N2_P2,'color','r','LineWidth', 1.5); hold on;
N2_P2 =squeeze(cond_1.data(:,chan_idx,:));
mean_N2_P2 = mean(N2_P2);

SEM=std(N2_P2)/sqrt(size(N2_P2,1));    %Standard Error
CI95=tinv([0.025 0.975],size(N2_P2,1)-1);                               %T-Score
ERP_CI95=bsxfun(@times, SEM,CI95(:));  %Confidence Intervals 
lower_b=mean_N2_P2+ERP_CI95(1,:);
upper_b=mean_N2_P2+ERP_CI95(2,:);
plot(my_time,mean_N2_P2,'color','b','LineWidth', 1.5); hold on;
plot(my_time,upper_b,'color','w')
plot(my_time,lower_b,'color','w'); 
x2=[my_time,fliplr(my_time)];
inbetween=[upper_b, fliplr(lower_b)];
fill(x2,inbetween,'b','FaceAlpha',0.2,'EdgeColor','none'); %set colour black, transparency 10%, and no edgecolour

N2_P2 =squeeze(cond_3.data(:,chan_idx,:));
mean_N2_P2 = mean(N2_P2);
SEM=std(N2_P2)/sqrt(size(N2_P2,1));    %Standard Error
CI95=tinv([0.025 0.975],size(N2_P2,1)-1);                               %T-Score
ERP_CI95=bsxfun(@times, SEM,CI95(:));  %Confidence Intervals 
lower_b=mean_N2_P2+ERP_CI95(1,:);
upper_b=mean_N2_P2+ERP_CI95(2,:);
plot(my_time,mean_N2_P2,'color','r','LineWidth', 1.5); hold on;
plot(my_time,upper_b,'color','w')
plot(my_time,lower_b,'color','w'); 
x2=[my_time,fliplr(my_time)];
inbetween=[upper_b, fliplr(lower_b)];
fill(x2,inbetween,'r','FaceAlpha',0.2,'EdgeColor','none'); %set colour black, transparency 10%, and no edgecolour
xlim(my_lims_N2P2);
ylim(my_xlims);
set(gca, 'YDir','reverse')
%xlabel('Times (ms)')
ylabel('\muV')
legend('Contact heat','Laser')

%% lowpass filtering

lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
         'PassbandFrequency',20,'PassbandRipple',0.2, ...
         'SampleRate',500);
data_in = cond_1.data;
n_subs = size(data_in,1);
n_chans = size(data_in,2);
filtered_data = zeros(size(data_in));

for i = 1:n_subs
    for k = 1:n_chans
        data = squeeze(data_in(i,k,:));
        filtered_data(i,k,:) = filtfilt(lpFilt,data);
    end
   chan_idx = 2;

end
cond_1.filt_data = filtered_data;
cond_1.filt_data(9,:,:) = [];

figure
chan_idx = 2;
data2plot = squeeze(cond_1.filt_data(:,chan_idx,:));
plot(data2plot')





data_in = cond_3.data;
n_subs = size(data_in,1);
n_chans = size(data_in,2);
filtered_data = zeros(size(data_in));

for i = 1:n_subs
    for k = 1:n_chans
        data = squeeze(data_in(i,k,:));
        filtered_data(i,k,:) = filtfilt(lpFilt,data);
    end

end
cond_3.filt_data = filtered_data;



%% making stim_parameters
laser_stim = [0;0;0;2.93;2.93;0;0];
laser_time = [-500 -200 0 1 5 6 1000];
% figure
% plot(laser_time,laser_stim);

contact_stim = [35;35;35;52;35;35];
contact_time = [-500 -250 0 245 670 1000];
% figure
% plot(contact_time,contact_stim);

my_lims_N2P2 = [-300 800];
fig = figure;
subplot(1,3,1)
yyaxis right
plot(contact_time,contact_stim);
ylabel('Temperature (C°)');
ylim([32 120]);
yticks([35 52])
yyaxis left
my_chans_in = {'c3','cz','c4','fz','pz'};
my_xlims = [-22 22];
chan_idx = 1;
N2_P2 =squeeze(cond_1.filt_data(:,chan_idx,:));
mean_N2_P2 = mean(N2_P2);
plot(my_time,mean_N2_P2,'color','k','LineWidth', 1.5); hold on;
SEM=std(N2_P2)/sqrt(size(N2_P2,1));    %Standard Error
CI95=tinv([0.025 0.975],size(N2_P2,1)-1);                               %T-Score
ERP_CI95=bsxfun(@times, SEM,CI95(:));  %Confidence Intervals 
lower_b=mean_N2_P2+ERP_CI95(1,:);
upper_b=mean_N2_P2+ERP_CI95(2,:);
plot(my_time,mean_N2_P2,'color','k','LineWidth', 1.5); hold on;
plot(my_time,upper_b,'color','w')
plot(my_time,lower_b,'color','w'); 
x2=[my_time,fliplr(my_time)];
inbetween=[upper_b, fliplr(lower_b)];
fill(x2,inbetween,'k','FaceAlpha',0.2,'EdgeColor','none'); %set colour black, transparency 10%, and no edgecolour
ylim(my_xlims);
ylabel('\muV')
set(gca, 'YDir','reverse')
xlim(my_lims_N2P2);
xlabel('Times (ms)')
xticks([-200 0 200 400 600])
set(gcf,'color','w');
subplot(1,3,2)
yyaxis right
plot(contact_time,contact_stim);
ylabel('Temperature (C°)');
ylim([32 120]);
yticks([35 52])
yyaxis left
my_chans_in = {'c3','cz','c4','fz','pz'};
my_xlims = [-22 22];
chan_idx = 2;
N2_P2 =squeeze(cond_1.filt_data(:,chan_idx,:));
mean_N2_P2 = mean(N2_P2);
plot(my_time,mean_N2_P2,'color','k','LineWidth', 1.5); hold on;
SEM=std(N2_P2)/sqrt(size(N2_P2,1));    %Standard Error
CI95=tinv([0.025 0.975],size(N2_P2,1)-1);                               %T-Score
ERP_CI95=bsxfun(@times, SEM,CI95(:));  %Confidence Intervals 
lower_b=mean_N2_P2+ERP_CI95(1,:);
upper_b=mean_N2_P2+ERP_CI95(2,:);
plot(my_time,mean_N2_P2,'color','k','LineWidth', 1.5); hold on;
plot(my_time,upper_b,'color','w')
plot(my_time,lower_b,'color','w'); 
x2=[my_time,fliplr(my_time)];
inbetween=[upper_b, fliplr(lower_b)];
fill(x2,inbetween,'k','FaceAlpha',0.2,'EdgeColor','none'); %set colour black, transparency 10%, and no edgecolour
ylim(my_xlims);
ylabel('\muV')
set(gca, 'YDir','reverse')
xlim(my_lims_N2P2);
xlabel('Times (ms)')
xticks([-200 0 200 400 600])
set(gcf,'color','w');
subplot(1,3,3)
yyaxis right
plot(contact_time,contact_stim);
ylabel('Temperature (C°)');
ylim([32 120]);
yticks([35 52])
yyaxis left
my_chans_in = {'c3','cz','c4','fz','pz'};
my_xlims = [-22 22];
chan_idx = 3;
N2_P2 =squeeze(cond_1.filt_data(:,chan_idx,:));
mean_N2_P2 = mean(N2_P2);
plot(my_time,mean_N2_P2,'color','k','LineWidth', 1.5); hold on;
SEM=std(N2_P2)/sqrt(size(N2_P2,1));    %Standard Error
CI95=tinv([0.025 0.975],size(N2_P2,1)-1);                               %T-Score
ERP_CI95=bsxfun(@times, SEM,CI95(:));  %Confidence Intervals 
lower_b=mean_N2_P2+ERP_CI95(1,:);
upper_b=mean_N2_P2+ERP_CI95(2,:);
plot(my_time,mean_N2_P2,'color','k','LineWidth', 1.5); hold on;
plot(my_time,upper_b,'color','w')
plot(my_time,lower_b,'color','w'); 
x2=[my_time,fliplr(my_time)];
inbetween=[upper_b, fliplr(lower_b)];
fill(x2,inbetween,'k','FaceAlpha',0.2,'EdgeColor','none'); %set colour black, transparency 10%, and no edgecolour
ylim(my_xlims);
ylabel('\muV')
set(gca, 'YDir','reverse')
xlim(my_lims_N2P2);
xlabel('Times (ms)')
xticks([-200 0 200 400 600])
set(gcf,'color','w');















fig = figure;
subplot(1,3,1)
yyaxis right
plot(laser_time,laser_stim);
ylabel('Intensity (J)');
ylim([-0.5 10]);
yticks([0 2])
yyaxis left
my_chans_in = {'c3','cz','c4','fz','pz'};
my_xlims = [-22 22];
chan_idx = 1;
N2_P2 =squeeze(cond_3.filt_data(:,chan_idx,:));
mean_N2_P2 = mean(N2_P2);
plot(my_time,mean_N2_P2,'color','k','LineWidth', 1.5); hold on;
SEM=std(N2_P2)/sqrt(size(N2_P2,1));    %Standard Error
CI95=tinv([0.025 0.975],size(N2_P2,1)-1);                               %T-Score
ERP_CI95=bsxfun(@times, SEM,CI95(:));  %Confidence Intervals 
lower_b=mean_N2_P2+ERP_CI95(1,:);
upper_b=mean_N2_P2+ERP_CI95(2,:);
plot(my_time,mean_N2_P2,'color','k','LineWidth', 1.5); hold on;
plot(my_time,upper_b,'color','w')
plot(my_time,lower_b,'color','w'); 
x2=[my_time,fliplr(my_time)];
inbetween=[upper_b, fliplr(lower_b)];
fill(x2,inbetween,'k','FaceAlpha',0.2,'EdgeColor','none'); %set colour black, transparency 10%, and no edgecolour
ylim(my_xlims);
ylabel('\muV')
set(gca, 'YDir','reverse')
xlim(my_lims_N2P2);
xlabel('Times (ms)')
xticks([-200 0 200 400 600])
set(gcf,'color','w');


subplot(1,3,2)
yyaxis right
plot(laser_time,laser_stim);
ylabel('Intensity (J)');
ylim([-0.5 10]);
yticks([0 2])
yyaxis left
my_chans_in = {'c3','cz','c4','fz','pz'};
my_xlims = [-22 22];
chan_idx = 2;
N2_P2 =squeeze(cond_3.filt_data(:,chan_idx,:));
mean_N2_P2 = mean(N2_P2);
plot(my_time,mean_N2_P2,'color','k','LineWidth', 1.5); hold on;
SEM=std(N2_P2)/sqrt(size(N2_P2,1));    %Standard Error
CI95=tinv([0.025 0.975],size(N2_P2,1)-1);                               %T-Score
ERP_CI95=bsxfun(@times, SEM,CI95(:));  %Confidence Intervals 
lower_b=mean_N2_P2+ERP_CI95(1,:);
upper_b=mean_N2_P2+ERP_CI95(2,:);
plot(my_time,mean_N2_P2,'color','k','LineWidth', 1.5); hold on;
plot(my_time,upper_b,'color','w')
plot(my_time,lower_b,'color','w'); 
x2=[my_time,fliplr(my_time)];
inbetween=[upper_b, fliplr(lower_b)];
fill(x2,inbetween,'k','FaceAlpha',0.2,'EdgeColor','none'); %set colour black, transparency 10%, and no edgecolour
ylim(my_xlims);
ylabel('\muV')
set(gca, 'YDir','reverse')
xlim(my_lims_N2P2);
xlabel('Times (ms)')
xticks([-200 0 200 400 600])
set(gcf,'color','w');

subplot(1,3,3)
yyaxis right
plot(laser_time,laser_stim);
ylabel('Intensity (J)');
ylim([-0.5 10]);
yticks([0 2])
yyaxis left
my_chans_in = {'c3','cz','c4','fz','pz'};
my_xlims = [-22 22];
chan_idx = 3;
N2_P2 =squeeze(cond_3.filt_data(:,chan_idx,:));
mean_N2_P2 = mean(N2_P2);
plot(my_time,mean_N2_P2,'color','k','LineWidth', 1.5); hold on;
SEM=std(N2_P2)/sqrt(size(N2_P2,1));    %Standard Error
CI95=tinv([0.025 0.975],size(N2_P2,1)-1);                               %T-Score
ERP_CI95=bsxfun(@times, SEM,CI95(:));  %Confidence Intervals 
lower_b=mean_N2_P2+ERP_CI95(1,:);
upper_b=mean_N2_P2+ERP_CI95(2,:);
plot(my_time,mean_N2_P2,'color','k','LineWidth', 1.5); hold on;
plot(my_time,upper_b,'color','w')
plot(my_time,lower_b,'color','w'); 
x2=[my_time,fliplr(my_time)];
inbetween=[upper_b, fliplr(lower_b)];
fill(x2,inbetween,'k','FaceAlpha',0.2,'EdgeColor','none'); %set colour black, transparency 10%, and no edgecolour
ylim(my_xlims);
ylabel('\muV')
set(gca, 'YDir','reverse')
xlim(my_lims_N2P2);
xlabel('Times (ms)')
xticks([-200 0 200 400 600])
set(gcf,'color','w');


%% Plotting


% Figure Cc - All conditions
my_clims = [-2.3 2.3];
x = cond_1.times;
y = cond_1.frex;
y(2) = 32;
y(4) = 37;

figure 
hold on
subplot(2,2,2)
data2plot = cond_1.Cc_mean;
contourf(x,y,data2plot,40,'linecolor','none');
title('Placebo Cheps');
caxis(my_clims)
cbar
% subplot(2,2,4)
% data2plot = cond_2.Cc_mean;
% contourf(x,y,data2plot,40,'linecolor','none');
% title('Cap Cheps');
% caxis(my_clims)
% cbar
subplot(2,2,1)
data2plot = cond_3.Cc_mean;
contourf(x,y,data2plot,40,'linecolor','none');
title('Placebo Leps');
caxis(my_clims)
cbar
% subplot(2,2,3)
% data2plot = cond_4.Cc_mean;
% contourf(x,y,data2plot,40,'linecolor','none');
% title('Cap Leps');
% caxis(my_clims)
% cbar



% Ci all conditions
my_clims = [-2 2];
figure
hold on
subplot(2,2,2)
data2plot = cond_1.Ci_mean;
contourf(x,y,data2plot,40,'linecolor','none');
title('Placebo Cheps');
caxis(my_clims)
cbar
subplot(2,2,4)
data2plot = cond_2.Ci_mean;
contourf(x,y,data2plot,40,'linecolor','none');
title('Cap Cheps');
caxis(my_clims)
cbar
subplot(2,2,1)
data2plot = cond_3.Ci_mean;
contourf(x,y,data2plot,40,'linecolor','none');
title('Placebo Leps');
caxis(my_clims)
cbar
subplot(2,2,3)
data2plot = cond_4.Ci_mean;
contourf(x,y,data2plot,40,'linecolor','none');
title('Cap Leps');
caxis(my_clims)
cbar



%% Plotting with missing subs

%Cond 2 - missing C08, C09
%Cond 1 - missing C02, C08

% subsetting the data
cond1_Ci = cond_1.Ci_gamma([1 3:end],:,:); % delete 7
cond2_Ci = cond_2.Ci_gamma([1 3:end],:,:); % delete 2
cond3_Ci = cond_3.Ci_gamma([1 4:7 9:end],:,:); % delete 2, 8 and 9
cond4_Ci = cond_4.Ci_gamma([1 3:7 10:end],:,:); % delete 2, 8 and 9

cond1_Cc = cond_1.Cc_gamma([1 3:end],:,:); 
cond2_Cc = cond_2.Cc_gamma([1 3:end],:,:);
cond3_Cc = cond_3.Cc_gamma([1 4:7 9:end],:,:);
cond4_Cc = cond_4.Cc_gamma([1 3:7 10:end],:,:);


% removing C03 also
cond1_Ci = cond_1.Ci_gamma([1:2 4:6 8:end],:,:); % delete 7
cond2_Ci = cond_2.Ci_gamma([1 4:end],:,:); % delete 2
cond3_Ci = cond_3.Ci_gamma([1 4:7 10:end],:,:); % delete 2, 8 and 9
cond4_Ci = cond_4.Ci_gamma([1 4:7 10:end],:,:); % delete 2, 8 and 9

cond1_Cc = cond_1.Cc_gamma([1:2 4:6 8:end],:,:); 
cond2_Cc = cond_2.Cc_gamma([1 4:end],:,:);
cond3_Cc = cond_3.Cc_gamma([1 4:7 10:end],:,:);
cond4_Cc = cond_4.Cc_gamma([1 4:7 10:end],:,:);




% computing new averages
avg_cond1_Ci = squeeze(mean(cond1_Ci));
 avg_cond2_Ci = squeeze(mean(cond2_Ci));
avg_cond3_Ci = squeeze(mean(cond3_Ci));
 avg_cond4_Ci = squeeze(mean(cond4_Ci));

avg_cond1_Cc = squeeze(mean(cond1_Cc));
avg_cond2_Cc = squeeze(mean(cond2_Cc));
avg_cond3_Cc = squeeze(mean(cond3_Cc));
 avg_cond4_Cc = squeeze(mean(cond4_Cc));

% Plotting again
% Figure Cc - All conditions
my_clims = [-2.3 2.3];
x = cond_1.times;
y = cond_1.frex;
y(2) = 32;
y(4) = 37;

f_h1 = figure;
hold on
subplot(1,2,2);
data2plot = avg_cond2_Cc;
contourf(x,y,data2plot,40,'linecolor','none');
title('Cc - Cap Cheps');
caxis(my_clims);
cbar;
% subplot(2,2,4);
% data2plot = avg_cond2_Cc;
% contourf(x,y,data2plot,40,'linecolor','none');
% title('Cc - Cap Cheps');
% caxis(my_clims);
% cbar;
subplot(1,2,1);
data2plot = avg_cond4_Cc;
contourf(x,y,data2plot,40,'linecolor','none');
title('Cc - Cap Leps');
caxis(my_clims);
cbar;
% subplot(2,2,3);
% data2plot = avg_cond4_Cc;
% contourf(x,y,data2plot,40,'linecolor','none');
% title('Cc - Cap Leps');
% caxis(my_clims);
% cbar;
% 
% 
% Ci
f_h2 = figure;
hold on
subplot(1,2,2)
data2plot = avg_cond1_Ci;
contourf(x,y,data2plot,40,'linecolor','none');
title('Ci - Placebo Cheps');
caxis(my_clims);
cbar;
% subplot(2,2,4);
% data2plot = avg_cond2_Ci;
% contourf(x,y,data2plot,40,'linecolor','none');
% title('Ci - Cap Cheps');
% caxis(my_clims);
% cbar;
subplot(1,2,1);
data2plot = avg_cond3_Ci;
contourf(x,y,data2plot,40,'linecolor','none');
title('Ci - Placebo Leps');
caxis(my_clims);
cbar;
% subplot(2,2,3);
% data2plot = avg_cond4_Ci;
% contourf(x,y,data2plot,40,'linecolor','none');
% title('Ci - Cap Leps');
% caxis(my_clims);
% cbar;


 saveas(f_h2,'Avg_Ci_no_S02_S03_S08_S09.jpg')
 saveas(f_h1,'Avg_Cc_no_S02_S03_S08_S09.jpg')

%% Plotting individual subject responses
S_ID = ['S01' ; 'S04' ; 'S05' ; 'S06' ; 'S07' ; 'S09' ; 'S10' ; 'S11';...
    'S12' ; 'S13' ; 'S14' ; 'S15' ; 'S16' ; 'S17' ; 'S18' ; 'S19' ; 'S20'];

 for i = 1:17
     % Getting subjects color limits
     all_data_Cc = [squeeze(cond1_Cc(i,:)) , squeeze(cond3_Cc(i,:))];
     all_data_Ci = [squeeze(cond1_Ci(i,:)) , squeeze(cond3_Ci(i,:))];
     my_max_Cc = max(abs(all_data_Cc));
     my_clims_Cc = [-my_max_Cc my_max_Cc];
     my_max_Ci = max(abs(all_data_Ci));
     my_clims_Ci = [-my_max_Ci my_max_Ci];
     
     Participant_ID = S_ID(i,:);   
     
     f_h = figure; % Cc
     my_clims = my_clims_Cc;
     x = cond_1.times;
     y = cond_1.frex;
     y(2) = 32;
     y(4) = 37;
     subplot(1,2,2)
     data2plot = squeeze(cond1_Cc(i,:,:));
     contourf(x,y,data2plot,40,'linecolor','none');
     title([Participant_ID ' Cc - Placebo Cheps']);
     caxis(my_clims);
     cbar;
%      subplot(2,2,4);
%      data2plot = squeeze(cond2_Cc(i,:,:));
%      contourf(x,y,data2plot,40,'linecolor','none');
%      title([Participant_ID ' Cc - Cap Cheps']);
%      caxis(my_clims);
%      cbar;
     subplot(1,2,1);
     data2plot = squeeze(cond3_Cc(i,:,:));
     contourf(x,y,data2plot,40,'linecolor','none');
     title([Participant_ID ' Cc - Placebo Leps']);
     caxis(my_clims);
     cbar;
%      subplot(2,2,3);
%      data2plot = squeeze(cond4_Cc(i,:,:));
%      contourf(x,y,data2plot,40,'linecolor','none');
%      title([Participant_ID ' Cc - Cap Leps']);
%      caxis(my_clims);
%      cbar;
     saveas(f_h, [Participant_ID '_Cc.jpg'])
     
     f_h = figure; % Ci
     my_clims = my_clims_Ci;
     x = cond_1.times;
     y = cond_1.frex;
     y(2) = 32;
     y(4) = 37;
     subplot(1,2,2)
     data2plot = squeeze(cond1_Ci(i,:,:));
     contourf(x,y,data2plot,40,'linecolor','none');
     title([Participant_ID ' Ci - Placebo Cheps']);
     caxis(my_clims);
     cbar;
%      subplot(2,2,4);
%      data2plot = squeeze(cond2_Ci(i,:,:));
%      contourf(x,y,data2plot,40,'linecolor','none');
%      title([Participant_ID ' Ci - Cap Cheps']);
%      caxis(my_clims);
%      cbar;
     subplot(1,2,1);
     data2plot = squeeze(cond3_Ci(i,:,:));
     contourf(x,y,data2plot,40,'linecolor','none');
     title([Participant_ID ' Ci - Placebo Leps']);
     caxis(my_clims);
     cbar;
%      subplot(2,2,3);
%      data2plot = squeeze(cond4_Ci(i,:,:));
%      contourf(x,y,data2plot,40,'linecolor','none');
%      title([Participant_ID ' Ci - Cap Leps']);
%      caxis(my_clims);
%      cbar;
     saveas(f_h, [Participant_ID '_Ci.jpg'])

     
 end
 


 
 %% Dem statistics prraap
 
 
 
 
 
%% Notes
%- C01_trial325J only has 10 stims
%- C02 trial 1 deleted epoch 11 for high emg noise levels
%- C02 trial 4 deleted epoch number one for high noise levels
%- C05 trial 1 delted epoch 21 for high deviation from baseline for a couple of channels
%- C05 trial 4 delted epoch 11 and 13 for high frontal activity resembling a blink
%- C08 trial 3 deleted epoch 16 due to large movement artifact
%- C08 trial 4 deleted epoch 18 due to resetting of aplifyer 
%- C12 trial 3 deleted epoch 1 due to high levels of drift across all channels
%- C13 trial 1 deleted epoch 15, 16, and 22 due to high levels of noise
%- C13 trial 3 delted epoch 2 due to drift of the central electrodes. there seems to be an odd
%    time-locked response on the central electrodes in this trial. check gamma data
%- C15 trial 2 deleted epoch 5 and 19 because they had motion artifacts
%- C16 trial 1 deleted epoch 3 due to high EMG response
%- C16 trial 2 deleted epoch 5 and 13 due to emg noise and drift
%- C17 condition 2 has a noisy right mastoid channel
%- C18 trial 1 deleted epoch 5 because it is noisy
%- C18 trial 2 - there might be some EMG burst that could be gamma but did not delete it from the
%data yet
%- C19 trial 1  deleted epoch 19 becuase of noise
%- C20 trial 1 deleted epoch 8 and 16 cause of high noise levels
%- C20 trial 2 deleted epoch 1 and 2 cause of high EMG artifact on right mastoid


%% Extracting low frequency

pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\2_cleaned_data'; % folder 

low_freq_detection(pathname)


%% Grouping of time-frequency data

% addpath C:\Users\Kip\Documents\eeglab2021.1 % add eeglab path
cd  Z:\30_Oscar_Ortiz\LEPs_CHEPs\scripts_functions
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\Raw
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\Raw\Participant_C06_Trial 1,3 & cond 2, 4'; % folder 
addpath  C:\Users\cassi\OneDrive\Desktop\Masters\EEG_and_Gamma\eeglab2021.1 % add eeglab path

% Condition 1
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\Low_freq\Cond_1';
cond_1 = grouping_data_lf(pathname);

% Condition 2
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\Low_freq\Cond_2';
cond_2 = grouping_data_lf(pathname);

% Condition 3
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\Low_freq\Cond_3';
cond_3 = grouping_data_lf(pathname);

% Condition 4
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\Low_freq\Cond_4';
cond_4 = grouping_data_lf(pathname);


%% Plotting with missing subs

%Cond 2 - missing C08, C09
%Cond 1 - missing C02, C08

% subsetting the data
cond1_Ci = cond_1.Ci_lf([1 4:end],:,:); % 2,3 these are deleting index
cond3_Ci = cond_3.Ci_lf([1 4:8 10:end],:,:); % 2,3,9

cond1_Cc = cond_1.Cc_lf([1 4:end],:,:); % 2,3
cond3_Cc = cond_3.Cc_lf([1 4:8 10:end],:,:); % 2,3,9



% computing new averages
avg_cond1_Ci = squeeze(mean(cond1_Ci));
avg_cond3_Ci = squeeze(mean(cond3_Ci));

avg_cond1_Cc = squeeze(mean(cond1_Cc));
avg_cond3_Cc = squeeze(mean(cond3_Cc));


% Plotting again
% Figure Cc - All conditions
my_clims = [-7 7];
x = cond_1.times;
y = cond_1.frex;
% y(2) = 32;
% y(4) = 37;

f_h1 = figure;
hold on
subplot(1,2,2);
data2plot = avg_cond1_Cc;
contourf(x,y,data2plot,40,'linecolor','none');
title('Cc - Placebo Cheps');
caxis(my_clims);
cbar;
subplot(1,2,1);
data2plot = avg_cond3_Cc;
contourf(x,y,data2plot,40,'linecolor','none');
title('Cc - Placebo Leps');
caxis(my_clims);
cbar;



% Ci
f_h2 = figure;
hold on
subplot(1,2,2)
data2plot = avg_cond1_Ci;
contourf(x,y,data2plot,40,'linecolor','none');
title('Ci - Placebo Cheps');
caxis(my_clims);
cbar;
subplot(1,2,1);
data2plot = avg_cond3_Ci;
contourf(x,y,data2plot,40,'linecolor','none');
title('Ci - Placebo Leps');
caxis(my_clims);
cbar;



 saveas(f_h2,'Avg_Ci_low_freq_no_S02_S03_S08_S09.jpg')
 saveas(f_h1,'Avg_Cc_low_freq_no_S02_S03_S08_S09.jpg')

%% Plotting individual subject responses
S_ID = ['S01' ; 'S03'; 'S04' ; 'S05' ; 'S06' ; 'S07' ; 'S10' ; 'S11';...
    'S12' ; 'S13' ; 'S14' ; 'S15' ; 'S16' ; 'S17' ; 'S18' ; 'S19' ; 'S20'];

 for i = 1:17
     % Getting subjects color limits
     all_data_Cc = [squeeze(cond1_Cc(i,:)),squeeze(cond3_Cc(i,:))];
     all_data_Ci = [squeeze(cond1_Ci(i,:)),squeeze(cond3_Ci(i,:))];
     my_max_Cc = max(abs(all_data_Cc));
     my_clims_Cc = [-my_max_Cc my_max_Cc];
     my_max_Ci = max(abs(all_data_Ci));
     my_clims_Ci = [-my_max_Ci my_max_Ci];
     
     Participant_ID = S_ID(i,:);   
     
     f_h = figure; % Cc
     my_clims = my_clims_Cc;
     x = cond_1.times;
     y = cond_1.frex;

     subplot(1,2,2)
     data2plot = squeeze(cond1_Cc(i,:,:));
     contourf(x,y,data2plot,40,'linecolor','none');
     title([Participant_ID ' Cc - Placebo Cheps']);
     caxis(my_clims);
     cbar;
     subplot(1,2,1);
     data2plot = squeeze(cond3_Cc(i,:,:));
     contourf(x,y,data2plot,40,'linecolor','none');
     title([Participant_ID ' Cc - Placebo Leps']);
     caxis(my_clims);
     cbar;
     saveas(f_h, [Participant_ID '_low_freq_Cc.jpg'])
     
     f_h = figure; % Ci
     my_clims = my_clims_Ci;
     x = cond_1.times;
     y = cond_1.frex;
     subplot(1,2,2)
     data2plot = squeeze(cond1_Ci(i,:,:));
     contourf(x,y,data2plot,40,'linecolor','none');
     title([Participant_ID ' Ci - Placebo Cheps']);
     caxis(my_clims);
     cbar;
     subplot(1,2,1);
     data2plot = squeeze(cond3_Ci(i,:,:));
     contourf(x,y,data2plot,40,'linecolor','none');
     title([Participant_ID ' Ci - Placebo Leps']);
     caxis(my_clims);
     cbar;
     saveas(f_h, [Participant_ID '_low_freq_Ci.jpg'])

     
 end
 
 %% Determinging regions ofthe gamma spectrum that are significantly greater than baseline 
 
 % getting baseline-uncorrected gamma spectra
 pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\2_cleaned_data'; % folder 
saving_folder = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\Gamma_not_baselined_corrected';
gamma_detection_no_baseline(pathname,saving_folder)


%% Grouping of time-frequency data

% Condition 1
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\Gamma_not_baselined_corrected\Cond1'
cond_1 = grouping_data_v2(pathname);

% Condition 3
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\Gamma_not_baselined_corrected\Cond3';
cond_3 = grouping_data_v2(pathname);

% Condition 2
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\Gamma_not_baselined_corrected\Cond2'
cond_2 = grouping_data_v2(pathname);

% Condition 4
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\3_sorted_data\gamma_all_extracted\Cond4';
cond_4 = grouping_data_v2(pathname);

%% Plotting all data 

%% Plotting

% re arranging S015 that had a stim on the other hand 
l_hand_cond1 = 14;
l_hand_cond3 = 15;

% Cond 1 
temp = cond_1.gamma_data;
temp(l_hand_cond1,chan1_IDX,:,:) = cond_1.gamma_data(l_hand_cond1,chan2_IDX,:,:);
temp(l_hand_cond1,chan2_IDX,:,:) = cond_1.gamma_data(l_hand_cond1,chan1_IDX,:,:);
cond_1_hand_corrected = cond_1;
cond_1_hand_corrected.gamma_data = temp;
cond_1_hand_corrected.gamma_mean = squeeze(mean(temp));

temp = cond_1.lf_data;
temp(l_hand_cond1,chan1_IDX,:,:) = cond_1.lf_data(l_hand_cond1,chan2_IDX,:,:);
temp(l_hand_cond1,chan2_IDX,:,:) = cond_1.lf_data(l_hand_cond1,chan1_IDX,:,:);
cond_1_hand_corrected.lf_data = temp;
cond_1_hand_corrected.lf_mean = squeeze(mean(temp));


% Cond 3
temp = cond_3.gamma_data;
temp(l_hand_cond3,chan1_IDX,:,:) = cond_3.gamma_data(l_hand_cond1,chan2_IDX,:,:);
temp(l_hand_cond3,chan2_IDX,:,:) = cond_3.gamma_data(l_hand_cond1,chan1_IDX,:,:);
cond_3_hand_corrected = cond_3;
cond_3_hand_corrected.gamma_data = temp;
cond_3_hand_corrected.gamma_mean = squeeze(mean(temp));

temp = cond_3.lf_data;
temp(l_hand_cond3,chan1_IDX,:,:) = cond_3.lf_data(l_hand_cond1,chan2_IDX,:,:);
temp(l_hand_cond3,chan2_IDX,:,:) = cond_3.lf_data(l_hand_cond1,chan1_IDX,:,:);
cond_3_hand_corrected.lf_data = temp;
cond_3_hand_corrected.lf_mean = squeeze(mean(temp));

cond_1 = cond_1_hand_corrected;
cond_3 = cond_3_hand_corrected;




chan_name1 = "C3";
chan_name2 = "C4";
chan_name3 = "Cz";
for k = 1:length(cond_1.chan_names_stored)
    C_chan = string(cond_1.chan_names_stored(k));
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

% Figure gamma_all conditions
my_clims = [-1.9 1.9];
x = cond_1.times;
y = cond_1.gamma_frex;

figure 
my_main_title = 'Gamma Response';
hold on
%Chan1
subplot(2,3,1)
data2plot = squeeze(cond_1.gamma_mean(chan1_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([chan_name1 'Placebo CHEPs']);
caxis(my_clims)
cbar
subplot(2,3,4)
data2plot = squeeze(cond_3.gamma_mean(chan1_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title(['Placebo LEPss']);
caxis(my_clims)
cbar
%Chan2
subplot(2,3,3)
data2plot = squeeze(cond_1.gamma_mean(chan2_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([chan_name2 'Placebo CHEPs']);
caxis(my_clims)
cbar
subplot(2,3,6)
data2plot = squeeze(cond_3.gamma_mean(chan2_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title(['Placebo LEPS']);
caxis(my_clims)
cbar
%Chan3
subplot(2,3,2)
data2plot = squeeze(cond_1.gamma_mean(chan3_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([my_main_title chan_name3 ' Placebo CHEPs']);
caxis(my_clims)
cbar
subplot(2,3,5)
data2plot = squeeze(cond_3.gamma_mean(chan3_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([' Placebo LEPs']);
caxis(my_clims)
cbar
%colormap(rbw_good_one_2);






% Figure lf all conditions
my_clims = [-8 8];
x = cond_1.times;
y = cond_1.lf_frex;

figure 
my_main_title = 'Low Frequency Response';
hold on
%Chan1
subplot(2,3,1)
data2plot = squeeze(cond_1.lf_mean(chan1_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([chan_name1 'Placebo CHEPs']);
caxis(my_clims)
cbar
subplot(2,3,4)
data2plot = squeeze(cond_3.lf_mean(chan1_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title(['Placebo LEPss']);
caxis(my_clims)
cbar
%Chan2
subplot(2,3,3)
data2plot = squeeze(cond_1.lf_mean(chan2_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([chan_name2 'Placebo CHEPs']);
caxis(my_clims)
cbar
subplot(2,3,6)
data2plot = squeeze(cond_3.lf_mean(chan2_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title(['Placebo LEPS']);
caxis(my_clims)
cbar
%Chan3
subplot(2,3,2)
data2plot = squeeze(cond_1.lf_mean(chan3_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([my_main_title chan_name3 ' Placebo CHEPs']);
caxis(my_clims)
cbar
subplot(2,3,5)
data2plot = squeeze(cond_3.lf_mean(chan3_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([' Placebo LEPs']);
caxis(my_clims)
cbar
%colormap(rbw_good_one_2);


%% Non parametric permutation using the function
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data
load('grouped_data_all_chans_extracted.mat')
struct_in = cond_1_hand_corrected;

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

input_mat = squeeze(struct_in.gamma_data(:,chan1_IDX,:,:));
n_permutes = 1000;
v_time = struct_in.times;
v_freq = struct_in.gamma_frex;
time_window = [-150 800];
plotting_input = 1;


%% condition 1
% C3
my_title = ['Cheps gamma ' char(chan_name1)];
input_mat = squeeze(struct_in.gamma_data(:,chan1_IDX,:,:));
[c1_zmap_c3, c1_zmapthresh_c3, c1_zmapthresh_plc_c3, c1_zmapthresh_clc_c3, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Cheps gamma ' char(chan_name2)];
input_mat = squeeze(struct_in.gamma_data(:,chan2_IDX,:,:));
[c1_zmap_c4, c1_zmapthresh_c4, c1_zmapthresh_plc_c4, c1_zmapthresh_clc_c4, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Cheps gamma ' char(chan_name3)];
input_mat = squeeze(struct_in.gamma_data(:,chan3_IDX,:,:));
[c1_zmap_cz, c1_zmapthresh_cz, c1_zmapthresh_plc_cz, c1_zmapthresh_clc_cz, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


%% condition 3
struct_in = cond_3_hand_corrected;

% C3
my_title = ['Leps gamma ' char(chan_name1)];
input_mat = squeeze(struct_in.gamma_data(:,chan1_IDX,:,:));
[c3_zmap_c3, c3_zmapthresh_c3, c3_zmapthresh_plc_c3, c3_zmapthresh_clc_c3, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Leps gamma ' char(chan_name2)];
input_mat = squeeze(struct_in.gamma_data(:,chan2_IDX,:,:));
[c3_zmap_c4, c3_zmapthresh_c4, c3_zmapthresh_plc_c4, c3_zmapthresh_clc_c4, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Leps gamma ' char(chan_name3)];
input_mat = squeeze(struct_in.gamma_data(:,chan3_IDX,:,:));
[c3_zmap_cz, c3_zmapthresh_cz, c3_zmapthresh_plc_cz, c3_zmapthresh_clc_cz, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%% Low frequency NPP
%% condition 1
% C3
my_title = ['Cheps gamma ' char(chan_name1)];
input_mat = squeeze(struct_in.gamma_data(:,chan1_IDX,:,:));
[c1_zmap_c3, c1_zmapthresh_c3, c1_zmapthresh_plc_c3, c1_zmapthresh_clc_c3, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Cheps gamma ' char(chan_name2)];
input_mat = squeeze(struct_in.gamma_data(:,chan2_IDX,:,:));
[c1_zmap_c4, c1_zmapthresh_c4, c1_zmapthresh_plc_c4, c1_zmapthresh_clc_c4, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Cheps gamma ' char(chan_name3)];
input_mat = squeeze(struct_in.gamma_data(:,chan3_IDX,:,:));
[c1_zmap_cz, c1_zmapthresh_cz, c1_zmapthresh_plc_cz, c1_zmapthresh_clc_cz, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


%% condition 3
struct_in = cond_3_hand_corrected;

% C3
my_title = ['Leps gamma ' char(chan_name1)];
input_mat = squeeze(struct_in.gamma_data(:,chan1_IDX,:,:));
[c3_zmap_c3, c3_zmapthresh_c3, c3_zmapthresh_plc_c3, c3_zmapthresh_clc_c3, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

% C4
my_title =  ['Leps gamma ' char(chan_name2)];
input_mat = squeeze(struct_in.gamma_data(:,chan2_IDX,:,:));
[c3_zmap_c4, c3_zmapthresh_c4, c3_zmapthresh_plc_c4, c3_zmapthresh_clc_c4, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%Cz
my_title =  ['Leps gamma ' char(chan_name3)];
input_mat = squeeze(struct_in.gamma_data(:,chan3_IDX,:,:));
[c3_zmap_cz, c3_zmapthresh_cz, c3_zmapthresh_plc_cz, c3_zmapthresh_clc_cz, tfv_time] = npp_stb(input_mat,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


%% plot actual responses
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data
gamma_in = load('nppc_gamma.mat');
lf_in = load('nppc_lf.mat');
x = gamma_in.x;
y = gamma_in.y;
clims_in = [-5 5];

figure
hold on
subplot(2,3,1)
title_in = ["C3","Cheps"];
zmap = gamma_in.c1_zmap_c3;
mask_in  = gamma_in.c1_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(2,3,2)
title_in = ["Cz","Cheps"];
zmap = gamma_in.c1_zmap_cz;
mask_in  = gamma_in.c1_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(2,3,3)
title_in = ["C4","Cheps"];
zmap = gamma_in.c1_zmap_c4;
mask_in  = gamma_in.c1_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)

subplot(2,3,4)
title_in = ["Leps"];
zmap = gamma_in.c3_zmap_c3;
mask_in  = gamma_in.c3_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(2,3,5)
zmap = gamma_in.c3_zmap_cz;
mask_in  = gamma_in.c3_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(2,3,6)
zmap = gamma_in.c3_zmap_c4;
mask_in  = gamma_in.c3_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
sgtitle('Gamma response with cluster correction')
 

%% Low frequency reponse
x = lf_in.x;
y = lf_in.y;
clims_in = [-10 10];

figure
hold on
subplot(2,3,1)
title_in = ["C3","Cheps"];
zmap = lf_in.c1_zmap_c3;
mask_in  = lf_in.c1_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(2,3,2)
title_in = ["Cz","Cheps"];
zmap = lf_in.c1_zmap_cz;
mask_in  = lf_in.c1_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(2,3,3)
title_in = ["C4","Cheps"];
zmap = lf_in.c1_zmap_c4;
mask_in  = lf_in.c1_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)

subplot(2,3,4)
title_in = ["Leps"];
zmap = lf_in.c3_zmap_c3;
mask_in  = lf_in.c3_zmapthresh_clc_c3;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
hold on
subplot(2,3,5)
zmap = lf_in.c3_zmap_cz;
mask_in  = lf_in.c3_zmapthresh_clc_cz;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
subplot(2,3,6)
zmap = lf_in.c3_zmap_c4;
mask_in  = lf_in.c3_zmapthresh_clc_c4;
plot_zmaps_with_outline(zmap,mask_in,x,y,clims_in,title_in)
sgtitle('Low frequency response with cluster correction')


%% do the point by point comparison between cond 1 and 3

addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data
load('grouped_data_all_chans_extracted.mat')
struct_in1 = cond_1_hand_corrected;
struct_in2 = cond_3_hand_corrected;
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
my_title = 'NPP t-test condition 3 - codition 1- C4 gamma';
[gamma_zmap_cz, gamma_zmapthresh_cz, gamma_zmapthresh_plc_cz, gamma_zmapthresh_clc_cz, tfv_time_cz, v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);




% lf
v_freq = struct_in1.lf_frex;
input_mat1 = squeeze(struct_in1.lf_data(:,chan1_IDX,:,:));
input_mat2 = squeeze(struct_in2.lf_data([1:7 9:end],chan1_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 1- C3 low frequency';
[lf_zmap_c3, lf_zmapthresh_c3, lf_zmapthresh_plc_c3, lf_zmapthresh_clc_c3, tfv_time_c3, lf_v_freq_c3] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

input_mat1 = squeeze(struct_in1.lf_data(:,chan2_IDX,:,:));
input_mat2 = squeeze(struct_in2.lf_data([1:7 9:end],chan2_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 1- C4 low frequency';
[lf_zmap_c4, lf_zmapthresh_c4, lf_zmapthresh_plc_c4, lf_zmapthresh_clc_c4, tfv_time_c4, lf_v_freq_c4] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);


input_mat1 = squeeze(struct_in1.lf_data(:,chan3_IDX,:,:));
input_mat2 = squeeze(struct_in2.lf_data([1:7 9:end],chan3_IDX,:,:));
my_title = 'NPP t-test condition 3 - codition 1- C4 low frequency';
[lf_zmap_cz, lf_zmapthresh_cz, lf_zmapthresh_plc_cz, lf_zmapthresh_clc_cz, tfv_time_cz, lf_v_freq_cz] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);

%% overlapping significant signals vs areas of significant differences
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data
gamma_in = load('nppc_gamma.mat');
lf_in = load('nppc_lf.mat');
my_results  = load('npp_results.mat');
load('grouped_data_all_chans_extracted.mat');
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data
load('grouped_data_all_chans_extracted.mat')
struct_in1 = cond_1_hand_corrected;
struct_in2 = cond_3_hand_corrected;
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

%C3 gamma
mask1 = gamma_in.c1_zmapthresh_clc_c3;
mask2 = gamma_in.c3_zmapthresh_clc_c3;
mask_difference = my_results.gamma_zmapthresh_clc_c3~=0;
mask_signal = mask1+mask2~=0;
mask_final = mask_signal+mask_difference==2;

% subset data based on the mask
% getting some parameters
time_window = [-150 800];
v_time = gamma_in.x;

time_s = dsearchn(v_time',time_window(1)); % start time for tf window
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
tfv_time  = v_time(time_s:time_e);
nTimepoints = numel(tfv_time);

chan_index = 1;
data_in1 = squeeze(cond_1_hand_corrected.gamma_data(:,chan_index,:,time_s:time_e));
data_out1 = zeros(size(data_in1));
data_in2 = squeeze(cond_3_hand_corrected.gamma_data([1:7 9:end],chan_index,:,time_s:time_e));
data_out2 = zeros(size(data_in2));
means_output = zeros(19,2);
for i = 1:19
    s_dat = squeeze(data_in1(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(1)
    data_out1(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_cc(i,1) = mean(temp);
    
    
    s_dat = squeeze(data_in2(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(1)
    data_out2(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_cc(i,2) = mean(temp);
end







%C4 gamma
mask1 = gamma_in.c1_zmapthresh_clc_c4;
mask2 = gamma_in.c3_zmapthresh_clc_c4;
mask_difference = my_results.gamma_zmapthresh_clc_c4~=0;
mask_signal = mask1+mask2~=0;
mask_final = mask_signal+mask_difference==2;

% subset data based on the mask
% getting some parameters
time_window = [-150 800];
v_time = gamma_in.x;

time_s = dsearchn(v_time',time_window(1)); % start time for tf window
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
tfv_time  = v_time(time_s:time_e);
nTimepoints = numel(tfv_time);

chan_index = 2;
data_in1 = squeeze(cond_1_hand_corrected.gamma_data(:,chan_index,:,time_s:time_e));
data_out1 = zeros(size(data_in1));
data_in2 = squeeze(cond_3_hand_corrected.gamma_data([1:7 9:end],chan_index,:,time_s:time_e));
data_out2 = zeros(size(data_in2));
means_output_ci = zeros(19,2);
for i = 1:19
    s_dat = squeeze(data_in1(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(1)
    data_out1(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_ci(i,1) = mean(temp);
    
    
    s_dat = squeeze(data_in2(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    pause(1)
    data_out2(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_ci(i,2) = mean(temp);
end


figure
subplot(1,2,1)
g = {'Cheps', 'Leps'};
boxplot(means_output_cc,g)
hold on
title('Cc');
    for j = 1:19
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

subplot(1,2,2)
g = {'Cheps', 'Leps'};
boxplot(means_output_ci,g)
hold on
title('Ci');
    for j = 1:19
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
sgtitle('Mean differences in gamma')


%% Round 2

%% Removing subjects that are not in the dataset
% Cond1 = no subject 8
% Cond3 = has all?

%removing s08 from cond 3
temp_copy = cond_3;
rmv_idx = 8; % 8 for subject number eight

g = temp_copy.gamma_data;
new_g = g([1:7,9:end],:,:,:);
new_mean_g = squeeze(mean(new_g));

l = temp_copy.lf_data;
new_l = l([1:7,9:end],:,:,:);
new_mean_l = squeeze(mean(new_l));

cond_3.gamma_data = new_g;
cond_3.gamma_mean = new_mean_g;
cond_3.lf_data = new_l;
cond_3.lf_mean = new_mean_l;


%% visualizing
% getting index
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

% Figure gamma_all conditions
my_clims = [-1.9 1.9];
x = cond_1.times;
y = cond_1.gamma_frex;

figure 
my_main_title = 'Gamma Response';
hold on
%Chan1
subplot(2,3,1)
data2plot = squeeze(cond_1.gamma_mean(chan1_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([chan_name1 'Placebo CHEPs']);
caxis(my_clims)
cbar
subplot(2,3,4)
data2plot = squeeze(cond_3.gamma_mean(chan1_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title(['Placebo LEPss']);
caxis(my_clims)
cbar
%Chan2
subplot(2,3,3)
data2plot = squeeze(cond_1.gamma_mean(chan2_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([chan_name2 'Placebo CHEPs']);
caxis(my_clims)
cbar
subplot(2,3,6)
data2plot = squeeze(cond_3.gamma_mean(chan2_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title(['Placebo LEPS']);
caxis(my_clims)
cbar
%Chan3
subplot(2,3,2)
data2plot = squeeze(cond_1.gamma_mean(chan3_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([my_main_title chan_name3 ' Placebo CHEPs']);
caxis(my_clims)
cbar
subplot(2,3,5)
data2plot = squeeze(cond_3.gamma_mean(chan3_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([' Placebo LEPs']);
caxis(my_clims)
cbar
%colormap(rbw_good_one_2);

%% figure 2 low frequency
% Figure gamma_all conditions
my_clims = [-7 7];
x = cond_1.times;
y = cond_1.lf_frex;

figure 
my_main_title = 'LF Response';
hold on
%Chan1
subplot(2,3,1)
data2plot = squeeze(cond_1.lf_mean(chan1_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([chan_name1 'Placebo CHEPs']);
caxis(my_clims)
cbar
subplot(2,3,4)
data2plot = squeeze(cond_3.lf_mean(chan1_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title(['Placebo LEPss']);
caxis(my_clims)
cbar
%Chan2
subplot(2,3,3)
data2plot = squeeze(cond_1.lf_mean(chan2_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([chan_name2 'Placebo CHEPs']);
caxis(my_clims)
cbar
subplot(2,3,6)
data2plot = squeeze(cond_3.lf_mean(chan2_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title(['Placebo LEPS']);
caxis(my_clims)
cbar
%Chan3
subplot(2,3,2)
data2plot = squeeze(cond_1.lf_mean(chan3_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([my_main_title chan_name3 ' Placebo CHEPs']);
caxis(my_clims)
cbar
subplot(2,3,5)
data2plot = squeeze(cond_3.lf_mean(chan3_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([' Placebo LEPs']);
caxis(my_clims)
cbar
%colormap(rbw_good_one_2);




% Figure gamma_all conditions
my_clims = [-1.9 1.9];
x = cond_1.times;
y = cond_1.lf_frex;


%% do individual subject responses
for i = 1:19
figure 
my_main_title = 'Gamma Response';
hold on
%Chan1
subplot(2,3,1)
data2plot = squeeze(cond_1.gamma_data(i,chan1_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([chan_name1 'Placebo CHEPs']);
%caxis(my_clims)
cbar
subplot(2,3,4)
data2plot = squeeze(cond_3.gamma_data(i,chan1_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title(['Placebo LEPss']);
%caxis(my_clims)
cbar
%Chan2
subplot(2,3,3)
data2plot = squeeze(cond_1.gamma_data(i,chan2_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([chan_name2 'Placebo CHEPs']);
%caxis(my_clims)
cbar
subplot(2,3,6)
data2plot = squeeze(cond_3.gamma_data(i,chan2_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title(['Placebo LEPS']);
%caxis(my_clims)
cbar
%Chan3
subplot(2,3,2)
data2plot = squeeze(cond_1.gamma_data(i,chan3_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([my_main_title chan_name3 ' Placebo CHEPs']);
%caxis(my_clims)
cbar
subplot(2,3,5)
data2plot = squeeze(cond_3.gamma_data(i,chan3_IDX,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title([' Placebo LEPs']);
caxis(my_clims)
cbar
%colormap(rbw_good_one_2);

end

%% Plotting each region
my_clims = [-1.9 1.9];
x = cond_1.times;
y = cond_1.gamma_frex;

chan_idx =1;
figure
data2plot = squeeze(cond_3.gamma_mean(chan_idx,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title(cond_1.chan_names_stored(chan_idx))
caxis(my_clims)
colormap jet
set(gcf,'color','w');
exportgraphics(gcf,'leps_c3.eps','ContentType','vector')

chan_idx =2;
figure
data2plot = squeeze(cond_3.gamma_mean(chan_idx,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title(cond_1.chan_names_stored(chan_idx))
caxis(my_clims)
colormap jet
set(gcf,'color','w');
exportgraphics(gcf,'leps_cz.eps','ContentType','vector')


chan_idx =3;
figure
data2plot = squeeze(cond_3.gamma_mean(chan_idx,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title(cond_1.chan_names_stored(chan_idx))
caxis(my_clims)
colormap jet
set(gcf,'color','w');
exportgraphics(gcf,'leps_c4.eps','ContentType','vector')


chan_idx =4;
figure
data2plot = squeeze(cond_3.gamma_mean(chan_idx,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title(cond_1.chan_names_stored(chan_idx))
caxis(my_clims)
colormap jet
set(gcf,'color','w');
exportgraphics(gcf,'leps_fz.eps','ContentType','vector')


chan_idx =5;
figure
data2plot = squeeze(cond_3.gamma_mean(chan_idx,:,:));
contourf(x,y,data2plot,40,'linecolor','none');
title(cond_1.chan_names_stored(chan_idx))
caxis(my_clims)
colormap jet
set(gcf,'color','w');
exportgraphics(gcf,'leps_pz.eps','ContentType','vector')

