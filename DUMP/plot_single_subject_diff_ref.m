function F1 = plot_single_subject_diff_ref(~,~,EEG_MAST,chans,S_ID,my_xlims,my_clims)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

chan_idx_v = zeros(length(chans),1);
% Find C3 idx~
 for k = 1:31
     for j = 1:length(chans)
         C_chan = string(EEG_MAST.chanlocs(k).labels);
         if C_chan == chans(j)
            chan_idx_v(j)=k;
         end
     end
 end

 %%

figure
 chan_id = 1;
% subplot(3,3,1)
% data2plot = squeeze(EEG_Fz.tf_data(chan_idx_v(chan_id),:,:));
% x = EEG_Fz.tf_time;
% y = EEG_Fz.tf_frex;
% contourf(x,y,data2plot,40,'linecolor','none');
% title([char(chans(chan_id)),'-Fz'])
% ylabel('Frequency (Hz)')
% xlim(my_xlims);
% clim(my_clims);
% colormap jet
% 
% subplot(3,3,4)
% data2plot = squeeze(EEG_AVG.tf_data(chan_idx_v(chan_id),:,:));
% x = EEG_Fz.tf_time;
% y = EEG_Fz.tf_frex;
% contourf(x,y,data2plot,40,'linecolor','none');
% title([char(chans(chan_id)),'-Avg'])
% ylabel('Frequency (Hz)')
% xlim(my_xlims);
% clim(my_clims);
% colormap jet


subplot(1,3,1)
data2plot = squeeze(EEG_MAST.tf_data(chan_idx_v(chan_id),:,:));
x = EEG_MAST.tf_time;
y = EEG_MAST.tf_frex;
contourf(x,y,data2plot,40,'linecolor','none');
title([char(chans(chan_id)),'-Mast'])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
xlim(my_xlims);
clim(my_clims);
colormap jet



chan_id = 2;
% subplot(3,3,2)
% data2plot = squeeze(EEG_Fz.tf_data(chan_idx_v(chan_id),:,:));
% x = EEG_Fz.tf_time;
% y = EEG_Fz.tf_frex;
% contourf(x,y,data2plot,40,'linecolor','none');
% title([char(chans(chan_id)),'-Fz'])
% xlim(my_xlims);
% clim(my_clims);
% colormap jet
% 
% subplot(3,3,5)
% data2plot = squeeze(EEG_AVG.tf_data(chan_idx_v(chan_id),:,:));
% x = EEG_Fz.tf_time;
% y = EEG_Fz.tf_frex;
% contourf(x,y,data2plot,40,'linecolor','none');
% title([char(chans(chan_id)),'-Avg'])
% xlim(my_xlims);
% clim(my_clims);
% colormap jet

subplot(1,3,2)
data2plot = squeeze(EEG_MAST.tf_data(chan_idx_v(chan_id),:,:));
x = EEG_MAST.tf_time;
y = EEG_MAST.tf_frex;
contourf(x,y,data2plot,40,'linecolor','none');
title([char(chans(chan_id)),'-Mast'])
xlabel('Time (ms)')
xlim(my_xlims);
clim(my_clims);
colormap jet


chan_id = 3;
% subplot(3,3,3)
% data2plot = squeeze(EEG_Fz.tf_data(chan_idx_v(chan_id),:,:));
% x = EEG_Fz.tf_time;
% y = EEG_Fz.tf_frex;
% contourf(x,y,data2plot,40,'linecolor','none');
% title([char(chans(chan_id)),'-Fz'])
% xlim(my_xlims);
% clim(my_clims);
% colormap jet
% %cbar
% 
% subplot(3,3,6)
% data2plot = squeeze(EEG_AVG.tf_data(chan_idx_v(chan_id),:,:));
% x = EEG_Fz.tf_time;
% y = EEG_Fz.tf_frex;
% contourf(x,y,data2plot,40,'linecolor','none');
% title([char(chans(chan_id)),'-Avg'])
% xlim(my_xlims);
% clim(my_clims);
% colormap jet
% %cbar

subplot(1,3,3)
data2plot = squeeze(EEG_MAST.tf_data(chan_idx_v(chan_id),:,:));
x = EEG_MAST.tf_time;
y = EEG_MAST.tf_frex;
contourf(x,y,data2plot,40,'linecolor','none');
title([char(chans(chan_id)),'-Mast'])
xlabel('Time (ms)')
xlim(my_xlims);
clim(my_clims);
colormap jet
cbar;
sgtitle(S_ID)

end