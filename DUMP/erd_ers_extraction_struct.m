function struct_out = erd_ers_extraction_struct(struct_in,time_windows_to_extract,length_window)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
struct_out = struct_in;

x = struct_in.tf_time;
x1 = struct_in.zmaps.time_v;
Y = struct_in.tf_frex;
n_chans = length(struct_in.chanlocs);
n_subs = size(struct_in.tf_data_subs,1);
n_time_w = length(time_windows_to_extract);
ERS_mean_v = zeros(n_subs,n_chans,n_time_w);
ERS_pp_v = ERS_mean_v;
ERS_w_mean_v = ERS_mean_v;
ERD_mean_v = ERS_mean_v;
ERD_pp_v = ERS_mean_v;
ERD_w_mean_v = ERS_mean_v;


for i = 1:n_chans
    data_in = squeeze(struct_in.tf_data_subs(:,i,:,:));
    mask_in = squeeze(struct_in.zmaps.zmapthresh_clc(i,:,:,:));
    for j = 1:n_subs
        data_s = squeeze(data_in(j,:,:));
        for n = 1:length(time_windows_to_extract)
            time_window = [time_windows_to_extract(n)-floor(length_window/2) time_windows_to_extract(n)+floor(length_window/2)];
            [ERS_mean_v(j,i,n),ERS_pp_v(j,i,n),ERS_w_mean_v(j,i,n),ERD_mean_v(j,i,n),ERD_pp_v(j,i,n),ERD_w_mean_v(j,i,n)] = erd_ers_extraction(mask_in,data_s,time_window,x,x1,Y);
       
        end
    end
end

struct_out.extraction_results.mean_ERS = ERS_mean_v;
struct_out.extraction_results.pp_ERS = ERS_pp_v;
struct_out.extraction_results.w_mean_ERS = ERS_w_mean_v;

struct_out.extraction_results.mean_ERD = ERD_mean_v;
struct_out.extraction_results.pp_ERD = ERD_pp_v;
struct_out.extraction_results.w_mean_ERD = ERD_w_mean_v;
struct_out.extraction_results.time_windows = time_windows_to_extract;


end

