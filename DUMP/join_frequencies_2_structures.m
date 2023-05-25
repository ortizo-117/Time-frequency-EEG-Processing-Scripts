function [struct_out1, struct_out2] = join_frequencies_2_structures(struct_in1,struct_in2)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
    struct_out1= struct_in1;
    struct_out2= struct_in2;
    gamma_1 = struct_in1.gamma_data;
    lf_1 = struct_in1.lf_data(:,:,1:end-1,:);
    gamma_m_1 = struct_in1.gamma_mean(:,:,:);

    gamma_2 = struct_in2.gamma_data;
    lf_2 = struct_in2.lf_data(:,:,1:end-1,:);
    gamma_m_2 = struct_in2.gamma_mean(:,:,:);
    
    gamma_m = [gamma_m_1 gamma_m_2];
    % getting normalized_outputs
    bound_gamma = max([max(max(gamma_m,[],2),[],3), abs(min(min(gamma_m,[],2),[],3))],[],2);
    inter_bound_gamma = ceil(bound_gamma);
    normalizing_fact_gamma  = max(inter_bound_gamma);
    norm_gamma_m_1 = gamma_m_1 ./ normalizing_fact_gamma;
    norm_gamma_m_2 = gamma_m_2 ./ normalizing_fact_gamma;
    
    % getting normalized_outputs
     lf_m_1 = struct_in1.lf_mean(:,1:end-1,:);
     lf_m_2 = struct_in2.lf_mean(:,1:end-1,:);
%     lf_m = [lf_m_1 lf_m_2]
    lim_factor = 4;
    lf_m_1_test =lf_m_1;
    lf_m_1_test(lf_m_1>lim_factor)= lim_factor;% change color limit scale
    lf_m_1_test(lf_m_1<-lim_factor)= -lim_factor;% change color limit scale

    lf_m_2_test =lf_m_2;
    lf_m_2_test(lf_m_2>lim_factor)= lim_factor;% change color limit scale
    lf_m_2_test(lf_m_2<-lim_factor)= -lim_factor;% change color limit scale
    lf_m = [lf_m_1_test lf_m_2_test];
    
    bound_lf = max([max(max(lf_m,[],2),[],3), abs(min(min(lf_m,[],2),[],3))],[],2);
    bounds_lf = [max(max(lf_m,[],2),[],3), abs(min(min(lf_m,[],2),[],3))];

    inter_bound_lf = ceil(bound_lf);
    normalizing_fact_lf  = max(inter_bound_lf);
    norm_lf_m_1 = lf_m_1 ./ normalizing_fact_lf;
    norm_lf_m_2 = lf_m_2 ./ normalizing_fact_lf;


    all_mean_1 = cat(2,lf_m_1,gamma_m_1);
    all_mean_norm_1 = cat(2,norm_lf_m_1,norm_gamma_m_1);
    all_fs = [struct_in1.lf_frex(1:end-1),struct_in1.gamma_frex];
    all_time = struct_in1.times;
    
    all_mean_2 = cat(2,lf_m_2,gamma_m_2);
    all_mean_norm_2 = cat(2,norm_lf_m_2,norm_gamma_m_2);


    %struct_out1.alldata = all_data;
    struct_out1.allmean_n = all_mean_norm_1;
    struct_out1.allmean = all_mean_1;
    struct_out1.allfrex = all_fs;
    struct_out1.all_time = all_time;
    struct_out1.gamma_bounds = bound_gamma;
    struct_out1.gamma_normalizing_fact = normalizing_fact_gamma;
    struct_out1.lf_bounds = bound_lf;
    struct_out1.lf_normalizing_fact = normalizing_fact_lf;

    %struct_out2.alldata = all_data_2;
    struct_out2.allmean_n = all_mean_norm_2;
    struct_out2.allmean = all_mean_2;
    struct_out2.allfrex = all_fs;
    struct_out2.all_time = all_time;
    struct_out2.gamma_bounds = bound_gamma;
    struct_out2.gamma_normalizing_fact = normalizing_fact_gamma;
    struct_out2.lf_bounds = bound_lf;
    struct_out2.lf_normalizing_fact = normalizing_fact_lf;

end