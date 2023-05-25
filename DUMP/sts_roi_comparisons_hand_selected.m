function [output,f1,f2] = sts_roi_comparisons_hand_selected(title_in,data_in1, data_in2,mask_sig1, mask_sig2,v_time,v_freq,v_time_results)
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


%    % get indices for 
%    [~,st_idx]=min(abs(v_time-time_w(1)));
%    [~,et_idx]=min(abs(v_time-time_w(2)));
%    [~,sf_idx]=min(abs(v_freq-freq_w(1)));
%    [~,ef_idx]=min(abs(v_freq-freq_w(2)));

   % subset data and masks into RoI as defined by time_w and freq_w
   D1 = fixed_data1(:,:,:);
   D2 = fixed_data2(:,:,:);
   M1 = mask_sig1(:,:);
   M2 = mask_sig2(:,:);
   
   D_vis = squeeze(mean(cat(1,D1,D2)));
   M_vis = logical(M1+M2);
   my_xlims = [0 800];
   clims_in = [-1.5 1.5];
   f1 = figure;
   tiledlayout(1,3)
   nexttile
   temp_title = 'Spectra of data 1';
   plot_og_dat_with_mask(squeeze(mean(D1)),M1,v_time_results,v_time_results,v_freq,clims_in,my_xlims,temp_title)
   nexttile
   temp_title = 'Combined spectra for the two electrodes';
   plot_og_dat_with_mask(D_vis,M_vis,v_time_results,v_time_results,v_freq,clims_in,my_xlims,temp_title)
   nexttile
   temp_title = 'Spectra of data 2';
   plot_og_dat_with_mask(squeeze(mean(D2)),M2,v_time_results,v_time_results,v_freq,clims_in,my_xlims,temp_title)
   %sgtitle(title_in);
   
   %% promt to determin ROIs and make calculations
   prompt='How many ROIs do you want to select?';
   n_RoI =  str2num((input(prompt,"s")));
   
   for i = 1:n_RoI
       name_roi = ['ROI_' num2str(i)];
       
       
       prompt=['What is the time (in ms) range (eg. 0 100) of RoI number ' num2str(i) '?'];
       t_range =  str2num((input(prompt,"s")));
       
       prompt=['What is the frequency (in Hz) range (eg. 30 100) of RoI number ' num2str(i) '?'];
       f_range =  str2num((input(prompt,"s")));
       

       
       [~,st_idx]=min(abs(v_time_results-t_range(1)));
       [~,et_idx]=min(abs(v_time_results-t_range(2)));
       [~,sf_idx]=min(abs(v_freq-f_range(1)));
       [~,ef_idx]=min(abs(v_freq-f_range(2)));
       
       subset1 = squeeze(mean(mean(D1(:,sf_idx:ef_idx,st_idx:et_idx),2),3));
       subset2 = squeeze(mean(mean(D2(:,sf_idx:ef_idx,st_idx:et_idx),2),3));
       data_out = [subset1 , subset2];
       
       
       figure
       hold on
       [h,p,ci,stats] = ttest(data_out(:,1),data_out(:,2));
       my_results.h = h;
       my_results.p = p;
       my_results.ci = ci;
       my_results.stats = stats;
       g = {'Data 1', 'Data 2'};
       boxplot([subset1,subset2],g);
       
       for j = 1:size(D1,1)
        y1 = data_out(j,:); % get first condition
        x1 = [1 2];
        lh = plot(x1,y1,'Color',[.05 .05 .05],...
            'LineWidth',0.1,...
            'Marker','.',...
            'MarkerEdgeColor','k',...
            'MarkerSize',12);
        lh.Color =[1,0,0,0.2];
       end
       title([name_roi ' difference'])
       
       eval(['output.' name_roi '.time_range = t_range;']);
       eval(['output.' name_roi '.frequency_range = f_range;']);
       eval(['output.' name_roi '.data = data_out;']);
       eval(['output.' name_roi '.results = my_results;']);
        
   end
   
   %% plotting new ROIs over the selection
   
   f2=figure;
   tiledlayout(1,2)
   nexttile
   temp_title = 'Spectra of data 1';
   plot_og_dat_with_mask(squeeze(mean(D1)),M1,v_time_results,v_time_results,v_freq,clims_in,my_xlims,temp_title)
   hold on
   for i = 1:n_RoI
       f_r = eval(['output.ROI_' num2str(i) '.frequency_range']);
       t_r = eval(['output.ROI_' num2str(i) '.time_range']);
%        %top_border
%        line(t_r,[f_r(2),f_r(2)] ,'LineWidth',4,'LineStyle','--','Color','r');
%        %bottom_border
%        line(t_r,[f_r(1),f_r(1)] ,'LineWidth',4,'LineStyle','--','Color','r');
%        %left_border
%        line([t_r(1) t_r(1)] ,f_r ,'LineWidth',4,'LineStyle','--','Color','r');
%        %right_border
%        line([t_r(2) t_r(2)],f_r ,'LineWidth',4,'LineStyle','--','Color','r');
        h = rectangle('Position', [t_r(1), f_r(1), t_r(2)-t_r(1), f_r(2)-f_r(1)], ...
                'Curvature', 0.2, ...
                'FaceColor', [192/255,192/255,192/255, 0.2], ...
                'EdgeColor', [0/255,0/255,0/255, 1]);
        
   
   end
   
   

   nexttile
   temp_title = 'Spectra of data 2';
   plot_og_dat_with_mask(squeeze(mean(D2)),M2,v_time_results,v_time_results,v_freq,clims_in,my_xlims,temp_title)
   hold on
   for i = 1:n_RoI
       f_r = eval(['output.ROI_' num2str(i) '.frequency_range']);
       t_r = eval(['output.ROI_' num2str(i) '.time_range']);
%        %top_border
%        line(t_r,[f_r(2),f_r(2)] ,'LineWidth',4,'LineStyle','--','Color','r');
%        %bottom_border
%        line(t_r,[f_r(1),f_r(1)] ,'LineWidth',4,'LineStyle','--','Color','r');
%        %left_border
%        line([t_r(1) t_r(1)] ,f_r ,'LineWidth',4,'LineStyle','--','Color','r');
%        %right_border
%        line([t_r(2) t_r(2)],f_r ,'LineWidth',4,'LineStyle','--','Color','r');
%        
       h = rectangle('Position', [t_r(1), f_r(1), t_r(2)-t_r(1), f_r(2)-f_r(1)], ...
                'Curvature', 0.2, ...
                'FaceColor', [192/255,192/255,192/255, 0.2], ...
                'EdgeColor', [0/255,0/255,0/255, 1]);
       
       
       
   end
   sgtitle(title_in);
   
   
   
   

   end
