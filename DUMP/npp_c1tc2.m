function [zmap, zmapthresh, zmapthresh_plc, zmapthresh_clc, tfv_time, v_freq] = npp_c1tc2(input_mat1,input_mat2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title)
%% Description:
% takes two time frequency matrices (dimensions = n_subs, frequencies, 
% times) as inputs and runs a non parametic permutation t-test. across both
% coditions. Inputs include input mat, number of permutations (n_permutes),
 % vectors for time and frequency(v_time and v_freq), limits in the time 
 % axis where you want to perform the analysis (time_window), and the
 % plotting_input variable, in which if is set to 1, it will plot the 
 % responses. The output consist of the z_map, the logical mask for the 
 % zmap uncorrected for multple comparisons (zmapthresh), the logical mask
 % with a pixel level correction for multiple comparisons (zmapthresh_plc),
 % and the logical threshold map with cluster correction for mutliple
 % comparisons (zmapthesh_clc). The last output is the new time vector
 % based on time_window.
 
 
 
%% Testing the function manually
% input_mat1 = squeeze(cond_1_hand_corrected.gamma_data(:,1,:,:)); 
% input_mat2 = squeeze(cond_3_hand_corrected.gamma_data(:,1,:,:));
% input_mat2 = input_mat2([1:7 9:end],:,:); % remove the extra sub
% 
% v_time = cond_1_hand_corrected.times;
% v_freq = cond_1_hand_corrected.gamma_frex;
% n_permutes = 100; % note: try to use 1000 or more permutations for real data
% time_window = [-150 800];
% plotting_input = 1;



%% checking number of subs in each matrix
temp1 = size(input_mat1);
n_subs1 = temp1(1);
temp2 = size(input_mat2);
n_subs2 = temp2(1);

if n_subs1 ~= n_subs2
    msg = ' Error, the number of subjects in the conditions is different. this function only works with equal number of subjects in each condition.';
    error(msg)
end



%% getting some parameters
time_s = dsearchn(v_time',time_window(1)); % start time for tf window
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
tfv_time  = v_time(time_s:time_e);
nTimepoints = numel(tfv_time);



% statistical significant values for voxel and cluster analysis
voxel_pval   = 0.05;
cluster_pval = 0.05;
mcc_voxel_pval = 0.05; % mcc = multiple comparisons correction
mcc_cluster_pval = 0.05;


%% get the data in a format for analyisis 
temp = squeeze(input_mat1(1,:,time_s:time_e));
tf_out_size = size(temp);
tf_vect_size = size(temp(:));
mat_1 = zeros(n_subs1,tf_vect_size(1));
mat_2 = zeros(n_subs2,tf_vect_size(1));

for i = 1:n_subs1
    % for mat 1
    curr = squeeze(input_mat1(i,:,time_s:time_e));
    curr_vect = curr(:);
    mat_1(i,:) = curr_vect';
  
    % for mat 2
    curr = squeeze(input_mat2(i,:,time_s:time_e));
    curr_vect = curr(:);
    mat_2(i,:) = curr_vect';
    
end

%% Running actual T test on data
t_mat = zeros(1,tf_vect_size(1));
for i = 1:tf_vect_size(1)
    points_1 = mat_1(:,i);
    points_2 = mat_2(:,i);
    [~,~,~,stats] = ttest(points_2,points_1);
    results = stats.tstat;
    t_mat(i) = results;
end

% reshaping into tf matrix
t_spec = reshape(t_mat,tf_out_size(1),tf_out_size(2));


%% running the premutations

% initialize null hypothesis matrices
permuted_maxvals = zeros(n_permutes,2,length(v_freq));
permuted_vals    = zeros(n_permutes,length(v_freq),nTimepoints);
max_clust_info   = zeros(n_permutes,1);
max_pixel_pvals = zeros(n_permutes,2);


% permutation loop 
cond1_mat = zeros(size(mat_1));
cond2_mat = zeros(size(mat_2));
curr_perm = zeros(1,tf_vect_size(1)); % initializing a variable to temporary store the results of current permutation before reshaping back into tf spectograph
x_wait = 0;
f = waitbar(x_wait,'Performing the permutations...','Name','npp_c1tc2',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);

for permi=1:n_permutes
     if getappdata(f,'canceling')
        delete(f)
        msg = 'Error, user terminated the function.';
        error(msg)
    end
    x_wait = permi/n_permutes;
    waitbar(x_wait,f,'performing the permutations...');
    
    conditions = randi([1, 2], n_subs1,1); %% create a vector that randomly assings a 1 or a 2 to each subject
    for j = 1:n_subs1
        if conditions(j) == 1 % if the subject is assigned condition 1
            cond1_sub = mat_1(j,:); % leave labling as is
            cond2_sub = mat_2(j,:);
        elseif conditions(j) == 2 % if the subject is assigned condition 2, switch the pre and post lable 
            cond1_sub = mat_2(j,:); %% switching lables
            cond2_sub = mat_1(j,:);
        end   
        cond1_mat(j,:) = cond1_sub; % creating the results matrix with randomly shuffled subjects for pre condition (cond1)
        cond2_mat(j,:) = cond2_sub; % creating the results matrix with randomly shuffled subjects for post condition (cond2)
    end
    
    % running point-by-point ttest for permutaiton
    for i = 1:tf_vect_size(1)
        cond1_points = cond1_mat(:,i);
        cond2_points = cond2_mat(:,i);
        [~,~,~,stats] = ttest(cond2_points,cond1_points);
        results = stats.tstat;
        curr_perm(i) = results;
    end
    
    % reshaping curr perm
    reshap_curr_perm = reshape (curr_perm,tf_out_size(1),tf_out_size(2));% reshaping
%     figure; - for troubleshooting
%     imagesc(reshap_curr_perm);
    permuted_vals(permi,:,:) = reshap_curr_perm; % saving permutation 
%    
    % saving max for correction for multiple comparisons
    max_pixel_pvals(permi,:) = [ min(reshap_curr_perm(:)) max(reshap_curr_perm(:)) ];

    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    % note that here, clusters were obtained by parametrically thresholding
    % the t-maps
    reshap_curr_perm(abs(reshap_curr_perm)<tinv(1-voxel_pval,n_subs1-1))=0;
    clustinfo = bwconncomp(reshap_curr_perm);
    max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % notes: cellfun is superfast, and the zero accounts for empty maps
end
delete(f)


%% getting zmaps
zmap = (t_spec-squeeze(mean(permuted_vals,1)))./squeeze(std(permuted_vals));
zmapthresh = zmap;
% appliying uncorrected threshold
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=false;
zmapthresh=logical(zmapthresh);


% apply pixel-level corrected threshold
lower_threshold = prctile(max_pixel_pvals(:,1),    mcc_voxel_pval*100/2);
upper_threshold = prctile(max_pixel_pvals(:,2),100-mcc_voxel_pval*100/2);
zmapthresh_plc = zmap;
mask1 = zmap>lower_threshold;
mask2 = zmap<upper_threshold;
mask = logical(mask1 + mask2-1);
zmapthresh_plc(mask) = 0;

%% apply cluster-level corrected threshold
zmapthresh_clc = zmap;
% apply uncorrected pixel-level threshold
zmapthresh_clc(abs(zmapthresh_clc)<norminv(1-voxel_pval))=false;
% find islands and remove those smaller than cluster size threshold
clustinfo = bwconncomp(zmapthresh_clc);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);
clust_threshold = prctile(max_clust_info,100-mcc_cluster_pval*100);

% identify clusters to remove
whichclusters2remove = find(clust_info<clust_threshold);

% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh_clc(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end

if plotting_input == 1
    figure
    my_lxlim = time_window(1)+10;
    my_hxlim = time_window(2);
    my_clims = [-5 5];
    subplot(2,2,1)
    contourf(tfv_time,v_freq,zmap,40,'linecolor','none')
    axis square
    set(gca,'clim',my_clims,'xlim',[my_lxlim my_hxlim])
    title('Unthresholded Z map')
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    colormap jet
    c = colorbar;
    c.Label.String = 'Z score';

    subplot(2,2,2)
    contourf(tfv_time,v_freq,zmap,40,'linecolor','none')
    hold on
    contour(tfv_time,v_freq,zmapthresh,1,'linecolor','k')
    c = colorbar;
    c.Label.String = 'Z score';
    axis square
    set(gca,'clim',my_clims,'xlim',[my_lxlim my_hxlim])
    title('Uncorrected threshold Z map')
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')


    subplot(2,2,3)
    contourf(tfv_time,v_freq,zmapthresh_plc,40,'linecolor','none')
    axis square
    set(gca,'clim',my_clims,'xlim',[my_lxlim my_hxlim])
    title('Pixel-corrected Z map')
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    c = colorbar;
    c.Label.String = 'Z score';

    subplot(2,2,4)
    contourf(tfv_time,v_freq,zmapthresh_clc,40,'linecolor','none')
    axis square
    set(gca,'clim',my_clims,'xlim',[my_lxlim my_hxlim])
    title('Cluster-corrected Z map')
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    c = colorbar;
    c.Label.String = 'Z score';
    %colormap(rbw_good_one_2)
    sgtitle(my_title)
end


end
