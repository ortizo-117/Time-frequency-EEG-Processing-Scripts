%% testing non parametric premutation testing
%% Author: Oscar Ortiz
 % Date: Nov/2021
 % Description: Runs a non-parametric permutated t test for a time freqency matrix that outputs
 % z-score maps that are corrected at a pixel or cluster level.
 % Inputs
 %    - condition_1 and condition 2: The two data matrices you want to compare (dimensions: number of subjects X
 %      number of frequencies X number of time bins).
 %    - v_time: vector delimiting the time stamps (in ms)
 %    - v_freq: vector delimiting the central frequencies extracted. 
 %    - n_permutes: number of permutations you want to run the data for. 
 %    - n_subs: number of subjects you are comparing


 %% Input section
condition_1 = cond1_Cc; 
condition_2 = cond3_Cc; 
v_time = cond_1.times;
v_freq = cond_1.frex;
v_freq(2) = 33.0; % might not need to do this for your own data
v_freq(4) = 37.0; % might not need to do this for your own data
n_permutes = 100; % note: try to use 1000 or more permutations for real data
n_subs = 17;

%% getting parameters for permutations
time_s = dsearchn(v_time',-150); % start time for tf window
time_e = dsearchn(v_time',800); % end time for tf window
tfv_time  = v_time(time_s:time_e);
nTimepoints = numel(tfv_time);



% statistical significant values for voxel and cluster analysis
voxel_pval   = 0.05;
cluster_pval = 0.05;
mcc_voxel_pval = 0.05; % mcc = multiple comparisons correction
mcc_cluster_pval = 0.05;



%% get the data in a format for analyisis 
temp = squeeze(condition_1(1,:,time_s:time_e));
tf_out_size = size(temp);

tf_vect_size = size(temp(:));
pre_mat = zeros(n_subs,tf_vect_size(1));
post_mat = zeros(n_subs,tf_vect_size(1));

for i = 1:n_subs
    % for pre
    curr = squeeze(condition_1(i,:,time_s:time_e));
    curr_vect = curr(:);
    pre_mat(i,:) = curr_vect';
  
    % for post
    curr = squeeze(condition_2(i,:,time_s:time_e));
    curr_vect = curr(:);
    post_mat(i,:) = curr_vect';
    
end

% temp = size(pre_mat(:));
% t_mat = zeros(temp(1),4);

%% Running actual T test on data
t_mat = zeros(1,tf_vect_size(1));
for i = 1:tf_vect_size(1)
    pre_points = pre_mat(:,i);
    post_points = post_mat(:,i);
    [~,~,~,stats] = ttest(post_points,pre_points);
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
cond1_mat = zeros(size(pre_mat));
cond2_mat = zeros(size(pre_mat));
curr_perm = zeros(1,tf_vect_size(1)); % initializing a variable to temporary store the results of current permutation before reshaping back into tf spectograph

for permi=1:n_permutes
    conditions = randi([1, 2], n_subs,1); %% create a vector that randomly assings a 1 or a 2 to each subject
    for j = 1:n_subs
        if conditions(j) == 1 % if the subject is assigned condition 1
            cond1_sub = pre_mat(j,:); % leave labling as is
            cond2_sub = post_mat(j,:);
        elseif conditions(j) == 2 % if the subject is assigned condition 2, switch the pre and post lable 
            cond1_sub = post_mat(j,:); %% switching lables
            cond2_sub = pre_mat(j,:);
        end   
        cond1_mat(j,:) = cond1_sub; % creating the results matrix with randomly shuffled subjects for pre condition (cond1)
        cond1_mat(j,:) = cond2_sub; % creating the results matrix with randomly shuffled subjects for post condition (cond2)
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
    reshap_curr_perm(abs(reshap_curr_perm)<tinv(1-voxel_pval,n_subs-1))=0;
    clustinfo = bwconncomp(reshap_curr_perm);
    max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % notes: cellfun is superfast, and the zero accounts for empty maps
end


%% getting zmaps
zmap = (t_spec-squeeze(mean(permuted_vals,1)))./squeeze(std(permuted_vals));

figure
my_lxlim = -140;
my_hxlim = 800;
subplot(2,2,1)
contourf(tfv_time,v_freq,zmap,40,'linecolor','none')
axis square
set(gca,'clim',[-5 5],'xlim',[my_lxlim my_hxlim])
title('Unthresholded Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
colormap jet
c = colorbar;
c.Label.String = 'Z score';

%% apply uncorrected threshold
subplot(2,2,2)
contourf(tfv_time,v_freq,zmap,40,'linecolor','none')
zmapthresh = zmap;
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=false;
zmapthresh=logical(zmapthresh);
hold on
contour(tfv_time,v_freq,zmapthresh,1,'linecolor','k')
c = colorbar;
c.Label.String = 'Z score';

axis square
set(gca,'clim',[-5 5],'xlim',[my_lxlim my_hxlim])
title('Uncorrected threshold Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')



%% apply pixel-level corrected threshold
lower_threshold = prctile(max_pixel_pvals(:,1),    mcc_voxel_pval*100/2);
upper_threshold = prctile(max_pixel_pvals(:,2),100-mcc_voxel_pval*100/2);

zmapthresh = zmap;
zmapthresh(zmapthresh>lower_threshold & zmapthresh<upper_threshold)=0; %% change to statistic dist not the z score?
subplot(2,2,3)
contourf(tfv_time,v_freq,zmapthresh,40,'linecolor','none')
axis square
set(gca,'clim',[-5 5],'xlim',[my_lxlim my_hxlim])
title('Pixel-corrected Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
c = colorbar;
c.Label.String = 'Z score';


%% apply cluster-level corrected threshold
zmapthresh = zmap;
% uncorrected pixel-level threshold
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
% find islands and remove those smaller than cluster size threshold
clustinfo = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);
clust_threshold = prctile(max_clust_info,100-mcc_cluster_pval*100);

% identify clusters to remove
whichclusters2remove = find(clust_info<clust_threshold);

% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end
subplot(2,2,4)
contourf(tfv_time,v_freq,zmapthresh,40,'linecolor','none')
axis square
set(gca,'clim',[-5 5],'xlim',[my_lxlim my_hxlim])
title('Cluster-corrected Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
c = colorbar;
c.Label.String = 'Z score';




%% getting zmaps - oscar correction
zmap = (t_spec-squeeze(mean(permuted_vals,1)))./squeeze(std(permuted_vals));

figure
my_lxlim = -140;
my_hxlim = 800;
subplot(2,2,1)
contourf(tfv_time,v_freq,zmap,40,'linecolor','none')
axis square
set(gca,'clim',[-5 5],'xlim',[my_lxlim my_hxlim])
title('Unthresholded Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
colormap jet
c = colorbar;
c.Label.String = 'Z score';

% apply uncorrected threshold
subplot(2,2,2)
contourf(tfv_time,v_freq,zmap,40,'linecolor','none')
zmapthresh = zmap;
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=false;
zmapthresh=logical(zmapthresh);
hold on
contour(tfv_time,v_freq,zmapthresh,1,'linecolor','k')
c = colorbar;
c.Label.String = 'Z score';

axis square
set(gca,'clim',[-5 5],'xlim',[my_lxlim my_hxlim])
title('Uncorrected threshold Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')



% apply pixel-level corrected threshold
lower_threshold = prctile(max_pixel_pvals(:,1),    mcc_voxel_pval*100/2);
upper_threshold = prctile(max_pixel_pvals(:,2),100-mcc_voxel_pval*100/2);
% lower_threshold = (lower_threshold-squeeze(mean(permuted_vals,1)))./squeeze(std(permuted_vals));
% upper_threshold = (upper_threshold-squeeze(mean(permuted_vals,1)))./squeeze(std(permuted_vals));

zmapthresh = zmap;
zmapthresh(t_spec>lower_threshold & t_spec<upper_threshold)=0; %% change to statistic dist not the z score?
subplot(2,2,3)
contourf(tfv_time,v_freq,zmapthresh,40,'linecolor','none')
axis square
set(gca,'clim',[-5 5],'xlim',[my_lxlim my_hxlim])
title('Pixel-corrected Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
c = colorbar;
c.Label.String = 'Z score';


% apply cluster-level corrected threshold
zmapthresh = zmap;
% uncorrected pixel-level threshold
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
% find islands and remove those smaller than cluster size threshold
clustinfo = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);
clust_threshold = prctile(max_clust_info,100-mcc_cluster_pval*100);

% identify clusters to remove
whichclusters2remove = find(clust_info<clust_threshold);

% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end
subplot(2,2,4)
contourf(tfv_time,v_freq,zmapthresh,40,'linecolor','none')
axis square
set(gca,'clim',[-5 5],'xlim',[my_lxlim my_hxlim])
title('Cluster-corrected Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
c = colorbar;
c.Label.String = 'Z score';

