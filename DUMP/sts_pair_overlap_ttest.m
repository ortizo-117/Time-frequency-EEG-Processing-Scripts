function [means_output,h,p,ci,stats] = sts_pair_overlap_ttest(data_in1, data_in2,mask_sig1, mask_sig2,mask_diff, time_window, v_time,v_timeresults, plotting_input)
%% Author: Oscar Ortiz
% Date: Feb 2022
% Description: this function takes two time frequency matrix (dimensions = 
% n_subs, frequencies, times), as well as the threshold maps exported from
% the npp_stb function for each of the signals (mask_sig1, mask_sig2,
% respectively) and the threshold maps experted from npp_c1tc2 (mask_diff)
% and extracts the subset of pixels in the t-f window defined in the time 
% domain with time_window parameter (also requires a time vector input)that
% are A)significantly greater than baseline and B) significantly different 
% from condition 1 vs conditions 2. Then it extracts the mean activity for 
% thissubset of pixels as a mean estimate of significantly different activity
% between conditions per subject and it performs a pairwise t-test between
% the two sets of subjects. The function exports a matrix containing two
% rows of these estimates per subject as well as all the statistc
% characteristics for a pairwise t-test (see ttest for details). Finally,
% if the plotting input is set to one, it will make a boxplot displaying
% your data and will add an asterisk between the conditions if there is a
% significantly difference between the two. 


%getting overlapping regions of interest
mask_signal = mask_sig1+mask_sig2~=0;
mask_final = mask_signal+mask_diff==2;

%% getting time parameters
time_s = dsearchn(v_time',time_window(1)); % start time for tf window
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
tfv_time  = v_time(time_s:time_e);

 % align time dimentions of data to match the dimensions of results
 [~,st_dat]=min(abs(v_time-min(v_time_results)));
 [~,et_dat]=min(abs(v_time-max(v_time_results)));
 fixed_data1 = data_in1(:,:,st_dat:et_dat);
 fixed_data2 = data_in2(:,:,st_dat:et_dat);


%% checking even number of subs
temp1 = size(data_in1);
temp2 = size(data_in2);
nsubs1 = temp1(1);
nsubs2 = temp2(1);

if nsubs1 ~= nsubs2
    msg = 'Number of subjects in each data structure is not the same. Please check the first dimension of each structure.';
    error(msg)
end


% creating output mats
data_out1 = zeros(temp1);
data_out2 = zeros(temp2);
means_output = zeros(nsubs1,2);

% subsetting each subject
for i = 1:nsubs1
    s_dat = squeeze(data_in1(i,:,:));
    s_dat(~mask_final)=NaN;
    %imagesc(s_dat)
    data_out1(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output(i,1) = mean(temp);
    
    
    s_dat = squeeze(data_in2(i,:,:));
    s_dat(~mask_final)=NaN;
    %imagesc(s_dat)
    data_out2(i,:,:) = s_dat;
    temp = s_dat(:);
    ind = isnan(temp);
    temp(ind) =  [];
    means_output(i,2) = mean(temp);
end


if plotting_input == 1
    figure
    hold on
    my_ylims = [floor(min(min(means_output))-.5) ceil(max(max(means_output))+.5)];
    g = {'Group 1', 'Group 2'};
    boxplot(means_output,g);
    [h,p,ci,stats] = ttest(means_output(:,1),means_output(:,2));
    for j = 1:nsubs1
        y1 = means_output(j,:); % get first condition
        x1 = [1 2];
        lh = plot(x1,y1,'Color',[.05 .05 .05],...
            'LineWidth',0.1,...
            'Marker','.',...
            'MarkerEdgeColor','k',...
            'MarkerSize',12);
        lh.Color =[1,0,0,0.2];
    end
    if p<0.05
        yt = get(gca, 'YTick');
        axis([xlim    0  ceil(max(yt)*1.2)])
        xt = get(gca, 'XTick');
        hold on
        plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
    end
    ylim(my_ylims);
end



end
