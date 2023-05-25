function [ERS_mean,ERS_pp,ERS_w_mean,ERD_mean,ERD_pp,ERD_w_mean] = erd_ers_extraction(mask_in,data_in,time_window,x,x1,Y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% mask_in
% data_in 
% time_window = [0 350];
% x = struct_in.tf_time;
% x1 = struct_in.zmaps.time_v;
%Y = struct_in.tf_frex;


mask_in_good = mask_in ~= 0;  
[~,start_mask] = min(abs(x1'-time_window(1)));
[~,end_mask] = min(abs(x1'-time_window(2)));
[~,start_dat] = min(abs(x'-time_window(1)));
[~,end_dat] = min(abs(x'-time_window(2)));

cut_mask = mask_in_good(:,start_mask:end_mask);
cut_dat = data_in(:,start_dat:end_dat);


X = x(start_dat:end_dat);
X1 = x1(start_mask:end_mask);


[Xq,Yq] = meshgrid(min(X):1:max(X),min(Y):max(Y));

[X1q,Y1q] = meshgrid(min(X1):1:max(X1),min(Y):max(Y));

% xq = min(X):1:max(X);
% yq = min(Y):max(Y);


cut_dat_resized= interp2(X,Y,cut_dat,Xq,Yq);
cut_mask_resized = interp2(X1,Y,double(cut_mask),X1q,Y1q)>=0.5;
cut_dat_cleaned = cut_mask_resized .* cut_dat_resized;
% 
% figure
% subplot(3,1,1)
% imagesc(cut_dat_resized)
% subplot(3,1,2)
% imagesc(cut_mask_resized)
% subplot(3,1,3)
% imagesc(cut_dat_cleaned)


%% ERS
ERS_mask = cut_dat_cleaned>0;
ERS_data = cut_dat_cleaned .*ERS_mask;
ERS_mean = mean(ERS_data(:));
ERS_pp = nnz(ERS_mask)/numel(ERS_mask);
ERS_w_mean = ERS_mean * ERS_pp;
%% ERD
ERD_mask = cut_dat_cleaned<0;
ERD_data = cut_dat_cleaned .*ERD_mask;
ERD_mean = mean(ERD_data(:));
ERD_pp = nnz(ERD_mask)/numel(ERD_mask);
ERD_w_mean = ERD_mean * ERD_pp;

end

