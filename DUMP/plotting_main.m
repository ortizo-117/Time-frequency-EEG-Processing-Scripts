%% Plotting all the main figures
%% Plotting the actual responses
%% plot actual responses
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\4_grouped_data
load('rbw_map.mat')
gamma_in = load('nppc_gamma.mat');
lf_in = load('nppc_lf.mat');


x = gamma_in.x;
y = gamma_in.y;
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
sgtitle('Gamma response with cluster correction')

% Making significant difference plots for gamma for C3, Cz, and C4
x = my_results.tfv_time;
y = my_results.v_freq_c3;
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
x = lf_in.x;
y = lf_in.y;
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
sgtitle('Low frequency response with cluster correction')

% difference
x = my_results.tfv_time;
y = my_results.lf_v_freq_c3;
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



















%% overlapping significant signals vs areas of significant differences-gamma
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
means_output_cc = zeros(19,2);
for i = 1:19
    s_dat = squeeze(data_in1(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    data_out1(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_cc(i,1) = mean(temp);
    
    
    s_dat = squeeze(data_in2(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
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
data_in1 = squeeze(cond_1_hand_corrected.gamma_data(:,chan2_IDX,:,time_s:time_e));
data_out1 = zeros(size(data_in1));
data_in2 = squeeze(cond_3_hand_corrected.gamma_data([1:7 9:end],chan2_IDX,:,time_s:time_e));
data_out2 = zeros(size(data_in2));
means_output_ci = zeros(19,2);
for i = 1:19
    s_dat = squeeze(data_in1(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    data_out1(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_ci(i,1) = mean(temp);
    
    
    s_dat = squeeze(data_in2(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    data_out2(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_ci(i,2) = mean(temp);
end
means_output_ci(9,:) = [];

figure
my_ylims = [-3.5 2];
subplot(1,2,1)
g = {'Cheps', 'Leps'};
boxplot(means_output_cc,g)
[h,p,ci,stats] = ttest(means_output_cc(:,1),means_output_cc(:,2))
if p<0.05
yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.2)])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
end

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
    ylim(my_ylims);
[h__cc,p_cc,ci_cc,stats_cc] = ttest(means_output_cc(:,1),means_output_cc(:,2))
ylabel('Mean Gamma ERSP (% difference from baseline)')
subplot(1,2,2)
g = {'Cheps', 'Leps'};
boxplot(means_output_ci,g)
[h,p,ci,stats] = ttest(means_output_ci(:,1),means_output_ci(:,2))
if p<0.05
yt = get(gca, 'YTick');
axis([xlim    0  ceil(max(yt)*1.2)])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
end
hold on
title('Ci');
    for j = 1:18
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
ylim(my_ylims);
sgtitle('Mean differences in gamma')



%% overlapping significant signals vs areas of significant differences
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
mask1 = lf_in.c1_zmapthresh_clc_c3;
mask2 = lf_in.c3_zmapthresh_clc_c3;
mask_difference_1 = my_results.lf_zmapthresh_clc_c3~=0;
mask_difference_2 = my_results.lf_zmapthresh_plc_c3~=0;
mask_difference = mask_difference_1+mask_difference_2~=0;
mask_signal = mask1+mask2~=0;
mask_final = mask_signal+mask_difference==2;

% subset data based on the mask
% getting some parameters
time_window = [-150 800];
v_time = lf_in.x;

time_s = dsearchn(v_time',time_window(1)); % start time for tf window
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
tfv_time  = v_time(time_s:time_e);
nTimepoints = numel(tfv_time);

chan_index = 1;
data_in1 = squeeze(cond_1_hand_corrected.lf_data(:,chan1_IDX,:,time_s:time_e));
data_out1 = zeros(size(data_in1));
data_in2 = squeeze(cond_3_hand_corrected.lf_data([1:7 9:end],chan1_IDX,:,time_s:time_e));
data_out2 = zeros(size(data_in2));
means_output = zeros(19,2);
for i = 1:19
    s_dat = squeeze(data_in1(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    data_out1(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_cc(i,1) = mean(temp);
    
    
    s_dat = squeeze(data_in2(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    data_out2(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_cc(i,2) = mean(temp);
end







%C4 gamma
mask1 = lf_in.c1_zmapthresh_clc_c4;
mask2 = lf_in.c3_zmapthresh_clc_c4;
mask_difference_1 = my_results.lf_zmapthresh_clc_c4~=0;
mask_difference_2 = my_results.lf_zmapthresh_plc_c4~=0;
mask_difference = mask_difference_1+mask_difference_2~=0;
mask_signal = mask1+mask2~=0;
mask_final = mask_signal+mask_difference==2;

% subset data based on the mask
% getting some parameters
time_window = [-150 800];
v_time = lf_in.x;

time_s = dsearchn(v_time',time_window(1)); % start time for tf window
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
tfv_time  = v_time(time_s:time_e);
nTimepoints = numel(tfv_time);

chan_index = 3;
data_in1 = squeeze(cond_1_hand_corrected.lf_data(:,chan2_IDX,:,time_s:time_e));
data_out1 = zeros(size(data_in1));
data_in2 = squeeze(cond_3_hand_corrected.lf_data([1:7 9:end],chan2_IDX,:,time_s:time_e));
data_out2 = zeros(size(data_in2));
means_output_ci = zeros(19,2);
for i = 1:19
    s_dat = squeeze(data_in1(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    data_out1(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_ci(i,1) = mean(temp);
    
    
    s_dat = squeeze(data_in2(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    data_out2(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_ci(i,2) = mean(temp);
end




%Cz gamma
mask1 = lf_in.c1_zmapthresh_clc_cz;
mask2 = lf_in.c3_zmapthresh_clc_cz;
mask_difference_1 = my_results.lf_zmapthresh_clc_cz~=0;
mask_difference_2 = my_results.lf_zmapthresh_plc_cz~=0;
mask_difference = mask_difference_1+mask_difference_2~=0;
mask_signal = mask1+mask2~=0;
mask_final = mask_signal+mask_difference==2;

% subset data based on the mask
% getting some parameters
time_window = [-150 800];
v_time = lf_in.x;

time_s = dsearchn(v_time',time_window(1)); % start time for tf window
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
tfv_time  = v_time(time_s:time_e);
nTimepoints = numel(tfv_time);

chan_index = 2;
data_in1 = squeeze(cond_1_hand_corrected.lf_data(:,chan3_IDX,:,time_s:time_e));
data_out1 = zeros(size(data_in1));
data_in2 = squeeze(cond_3_hand_corrected.lf_data([1:7 9:end],chan3_IDX,:,time_s:time_e));
data_out2 = zeros(size(data_in2));
means_output_cz = zeros(19,2);
for i = 1:19
    s_dat = squeeze(data_in1(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    data_out1(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_cz(i,1) = mean(temp);
    
    
    s_dat = squeeze(data_in2(i,:,:));
    s_dat(~mask_final)=NaN;
    imagesc(s_dat)
    data_out2(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output_cz(i,2) = mean(temp);
end


figure
subplot(1,3,1)
my_lims = [-1.5 4];
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
    ylim(my_lims);

[h,p,ci,stats] = ttest(means_output_cc(:,1),means_output_cc(:,2))

subplot(1,3,2)
g = {'Cheps', 'Leps'};
boxplot(means_output_cz,g)
hold on
title('Cz');
    for j = 1:19
        y1 = means_output_cz(j,:); % get first condition
        x1 = [1 2];
        lh = plot(x1,y1,'Color',[.05 .05 .05],...
            'LineWidth',0.1,...
            'Marker','.',...
            'MarkerEdgeColor','k',...
            'MarkerSize',12);
        lh.Color =[1,0,0,0.2];
    end
    ylim(my_lims);

[h,p,ci,stats] = ttest(means_output_cz(:,1),means_output_cz(:,2))

subplot(1,3,3)
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
    ylim(my_lims);

[h,p,ci,stats] = ttest(means_output_ci(:,1),means_output_ci(:,2))
sgtitle('Mean differences in lf')

