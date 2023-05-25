function [means_output,npix_output,w_means_output,mean_results,n_pix_results,w_mean_results] = sts_roi_comparisons(title_in,data_in1, data_in2,mask_sig1, mask_sig2,norp, time_w,freq_w,v_time,v_freq,v_time_results,plotting_input)
% Author: Oscar Ortiz
% Date: 07/13/2022
% Purpose: To extract values specified in a certain window in a tf spectra
   
% % Checking for arguments
%    arguments
%         data_in1 (:,:,:) {mustBeNumeric} 
%         data_in2 (:,:,:) {mustBeNumeric} 
%         mask_sig1 (:,:) logical
%         mask_sig2 (:,:) logical
%         norp (:,:) char {mustBeMember(norp,{'positive','negative'})} = 'positive'
%         time_w (:,:) {mustBeNumeric} 
%         freq_w (:,:) {mustBeNumeric} 
%         v_time (:,:) {mustBeNumeric} 
%         v_freq (:,:) {mustBeNumeric} 
%    end


   %% subset the data into the window of interest

   % align time dimentions of data to match the dimensions of results
   [~,st_dat]=min(abs(v_time-min(v_time_results)));
   [~,et_dat]=min(abs(v_time-max(v_time_results)));
   fixed_data1 = data_in1(:,:,st_dat:et_dat);
   fixed_data2 = data_in2(:,:,st_dat:et_dat);


   % get indices for 
   [~,st_idx]=min(abs(v_time-time_w(1)));
   [~,et_idx]=min(abs(v_time-time_w(2)));
   [~,sf_idx]=min(abs(v_freq-freq_w(1)));
   [~,ef_idx]=min(abs(v_freq-freq_w(2)));

   % subset data and masks into RoI as defined by time_w and freq_w
   D1 = fixed_data1(:,sf_idx:ef_idx,st_idx:et_idx);
   D2 = fixed_data2(:,sf_idx:ef_idx,st_idx:et_idx);
   M1 = mask_sig1(sf_idx:ef_idx,st_idx:et_idx);
   M2 = mask_sig2(sf_idx:ef_idx,st_idx:et_idx);

   means_output = zeros(size(D1,1),2);
   npix_output = means_output;

   for i = 1:size(D1,1)
       %% check which portions within that window are significant and in the directionality of interest
       % check if you want positive or negative value;
       s_data1 = squeeze(D1(i,:,:));
       s_data2 = squeeze(D2(i,:,:));

       size_mat1 = length(s_data1(:));
       size_mat2 = length(s_data2(:));

       s_data1(~M1) = NaN;
       s_data2(~M2) = NaN;

       if strcmp(norp,'positive')
            s_data1(s_data1<0) = NaN;
            s_data2(s_data2<0) = NaN;
       elseif strcmp(norp,'negative')
            s_data1(s_data1>0) = NaN;
            s_data2(s_data2>0) = NaN;
       end


       % getting the mean per subject 
       s_1 = s_data1(:);
       s_1(isnan(s_1)) = [];
       s_1mean = mean(s_1);
       s_1pix = length(s_1);


       s_2 = s_data2(:);
       s_2(isnan(s_2)) = [];
       s_2mean = mean(s_2);
       s_2pix = length(s_2);



       means_output(i,1) = s_1mean;
       means_output(i,2) = s_2mean;
       npix_output(i,1) = s_1pix/size_mat1;
       npix_output(i,2) = s_2pix/size_mat2;

   end

   w_means_output = npix_output .* means_output;


   if plotting_input == 1
    figure
    hold on
    my_ylims = [floor(min(min(means_output))-.5) ceil(max(max(means_output))+.5)];
    g = {'Group 1', 'Group 2'};
    boxplot(means_output,g);
    [h,p,ci,stats] = ttest(means_output(:,1),means_output(:,2));
    for j = 1:size(D1,1)
        y1 = means_output(j,:); % get first condition
        x1 = [1 2];
        lh = plot(x1,y1,'Color',[.05 .05 .05],...
            'LineWidth',0.1,...
            'Marker','.',...
            'MarkerEdgeColor','k',...
            'MarkerSize',12);
        lh.Color =[1,0,0,0.2];
    end
%     if p<0.05
%         yt = get(gca, 'YTick');
%         axis([xlim    floor(min(yt)*1.2)  ceil(max(yt)*1.2)])
%         xt = get(gca, 'XTick');
%         hold on
%         plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
%     end
   % ylim(my_ylims);
    title([title_in 'mean ' norp ' activity'])
    mean_results.h = h;
    mean_results.p = p;
    mean_results.ci = ci;
    mean_results.stats = stats;

    % N pix
    figure
    hold on
    my_ylims = [floor(min(min(npix_output))-.5) ceil(max(max(npix_output))+.5)];
    g = {'Group 1', 'Group 2'};
    boxplot(npix_output,g);
    [h,p,ci,stats] = ttest(npix_output(:,1),npix_output(:,2));
    for j = 1:size(D1,1)
        y1 = npix_output(j,:); % get first sub
        x1 = [1 2];
        lh = plot(x1,y1,'Color',[.05 .05 .05],...
            'LineWidth',0.1,...
            'Marker','.',...
            'MarkerEdgeColor','k',...
            'MarkerSize',12);
        lh.Color =[1,0,0,0.2];
    end
%     if p<0.05
%         yt = get(gca, 'YTick');
%         axis([xlim    floor(min(yt)*1.2)  ceil(max(yt)*1.2)])
%         xt = get(gca, 'XTick');
%         hold on
%         plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
%     end
    %ylim(my_ylims);
    title([title_in 'percentage of ' norp ' pixels in RoI'])
    n_pix_results.h = h;
    n_pix_results.p = p;
    n_pix_results.ci = ci;
    n_pix_results.stats = stats;

    % weighted mean
    figure
    hold on
    my_ylims = [floor(min(min(w_means_output))-.5) ceil(max(max(w_means_output))+.5)];
    g = {'Group 1', 'Group 2'};
    boxplot(w_means_output,g);
    [h,p,ci,stats] = ttest(w_means_output(:,1),w_means_output(:,2));
    for j = 1:size(D1,1)
        y1 = w_means_output(j,:); % get first condition
        x1 = [1 2];
        lh = plot(x1,y1,'Color',[.05 .05 .05],...
            'LineWidth',0.1,...
            'Marker','.',...
            'MarkerEdgeColor','k',...
            'MarkerSize',12);
        lh.Color =[1,0,0,0.2];
    end
%     if p<0.05
%         yt = get(gca, 'YTick');
%         axis([xlim    floor(min(yt)*1.2)  ceil(max(yt)*1.2)])
%         xt = get(gca, 'XTick');
%         hold on
%         plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
%     end
    %ylim(my_ylims);
    title([ title_in 'weighted mean ' norp ' activity (mean X percent of pixels)'])
    w_mean_results.h = h;
    w_mean_results.p = p;
    w_mean_results.ci = ci;
    w_mean_results.stats = stats;

   end




end

