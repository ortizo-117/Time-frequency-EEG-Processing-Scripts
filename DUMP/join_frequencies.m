function struct_out = join_frequencies(struct_in)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
    struct_out= struct_in;
    gamma = struct_in.gamma_data;
    lf = struct_in.lf_data(:,:,1:end-1,:);
    all_data = cat(3,lf,gamma);
    gamma_m = struct_in.gamma_mean(:,:,:);
    % getting normalized_outputs
    bound_gamma = max([max(max(gamma_m,[],2),[],3), abs(min(min(gamma_m,[],2),[],3))],[],2);
    norm_gamma_m = gamma_m ./ bound_gamma;
    
    % getting normalized_outputs
    lf_m = struct_in.lf_mean(:,1:end-1,:);
    bound_lf = max([max(max(lf_m,[],2),[],3), abs(min(min(lf_m,[],2),[],3))],[],2);
    norm_lf_m = lf_m ./ bound_lf;


    all_mean = cat(2,lf_m,gamma_m);
    all_mean_norm = cat(2,norm_lf_m,norm_gamma_m);
    all_fs = [struct_in.lf_frex(1:end-1),struct_in.gamma_frex];
    all_time = struct_in.times;

    struct_out.alldata = all_data;
    struct_out.allmean_n = all_mean_norm;
    struct_out.allmean = all_mean;
    struct_out.allfrex = all_fs;
    struct_out.all_time = all_time;
    struct_out.gamma_bounds = bound_gamma;
    struct_out.lf_bounds = bound_lf;
end