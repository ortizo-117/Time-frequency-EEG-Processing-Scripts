%% testing non parametric premutation testing baseline vs post
%% Author: Oscar Ortiz
 % Date: Jan/2022
 % Description: Runs a non-parametric permutated t test for a time freqency matrix that outputs
 % z-score maps that are corrected at a pixel or cluster level.It compares
 % the activity before the onset of stimulus against the rest of the time
 % frequency window
 % Inputs
 %    - input_mat = (dimensions: number of subjects X number of frequencies X number of time bins).
 %    - v_time: vector delimiting the time stamps (in ms)
 %    - v_freq: vector delimiting the central frequencies extracted. 
 %    - n_permutes: number of permutations you want to run the data for. 
 %    - n_subs: number of subjects you are comparing


 %% Input section
input_mat = cond3_Cc; % inpuit mat
avg_mat = squeeze(mean(input_mat));
v_time = cond_1.times;
v_freq = cond_1.frex;
v_freq(2) = 33.0; % might not need to do this for your own data
v_freq(4) = 37.0; % might not need to do this for your own data
n_permutes = 5000; % note: try to use 1000 or more permutations for real data
n_subs = 17;

%% getting parameters for permutations
time_s = dsearchn(v_time',-150); % start time for tf window
time_0 = dsearchn(v_time',0);
time_e = dsearchn(v_time',800); % end time for tf window
tfv_time  = v_time(time_s:time_e);
nTimepoints = numel(tfv_time);





% statistical significant values for voxel and cluster analysis
voxel_pval   = 0.05;
cluster_pval = 0.05;
mcc_voxel_pval = 0.05; % mcc = multiple comparisons correction
mcc_cluster_pval = 0.05;



%% get the data in a vector format for analyisis 
temp =avg_mat(:,time_s:time_e);
tf_out_size = size(temp);

tf_vect_size = size(temp(:));
% pre_mat = zeros(n_subs,tf_vect_size(1));
% post_mat = zeros(n_subs,tf_vect_size(1));

% for i = 1:n_subs
%     % for pre
%     curr = squeeze(condition_1(i,:,time_s:time_e));
%     curr_vect = curr(:);
%     pre_mat(i,:) = curr_vect';
%   
%     % for post
%     curr = squeeze(condition_2(i,:,time_s:time_e));
%     curr_vect = curr(:);
%     post_mat(i,:) = curr_vect';
%     
% end

% temp = size(pre_mat(:));
% t_mat = zeros(temp(1),4);

%% Running actual T test on data

%% need to run it per frequency bin

mat_size = size(avg_mat(:,time_s:time_e));
t_mat = zeros(mat_size);

for i = 1:mat_size(1) %run per frequency bin
    prestim_estimate = mean(input_mat(:,i,time_s:time_0),3);
    poststim_points = squeeze(input_mat(:,i,time_s:time_e));
    prestim_points = repmat(prestim_estimate, 1,length(poststim_points));
    [~,~,~,stats] = ttest(poststim_points,prestim_points);
    results = stats.tstat;
    t_mat(i,:) = results;
end
% 
% % reshaping into tf matrix % DONT tHINK i   NEED IT
% t_spec = reshape(t_mat,tf_out_size(1),tf_out_size(2));

%% running the premutations

% initialize null hypothesis matrices
permuted_maxvals = zeros(n_permutes,2,length(v_freq));
max_clust_info   = zeros(n_permutes,1);
max_pixel_pvals = zeros(n_permutes,2);
%% permutation loop 
curr_perm = zeros(mat_size); % initializing a variable to temporary store the results of current permutation before reshaping back into tf spectograph
permuted_vals = zeros([n_permutes, mat_size]);
for permi=1:n_permutes
    for i = 1:mat_size(1) % per frequency bin 
        prestim_estimate = squeeze(input_mat(:,i,time_s:time_0));
        poststim_points = squeeze(input_mat(:,i,time_s:time_e));
        pre_stim_estimate_r = zeros([n_subs,1]);
        for m = 1:n_subs
            r_index = randi([1,55]);
            pre_stim_estimate_r(m) = prestim_estimate(m,r_index);
        end
        prestim_points = repmat(pre_stim_estimate_r, 1,length(poststim_points));
        shuffled_data_1 = zeros(size(poststim_points));
        shuffled_data_2 = zeros(size(poststim_points));
        for j = 1:mat_size(2) % per timepoint
            c_data = poststim_points(:,j);
            c_data_out_1 = zeros(size(c_data));
            c_baseline = pre_stim_estimate_r;
            c_data_out_2 = zeros(size(c_data));
            conditions = randi([1, 2], n_subs,1); %% create a vector that randomly assings a 1 or a 2 to each subject
            for k = 1:n_subs % switching labels per subject
                if conditions(k) == 1 % if the subject is assigned condition 1
                    c_data_out_1(k) = c_data(k);
                    c_data_out_2(k) = c_baseline(k);
                elseif conditions(k) == 2
                    c_data_out_1(k) = c_baseline(k);
                    c_data_out_2(k) = c_data(k);
                end
            end
            shuffled_data_1(:,j) = c_data_out_1;
            shuffled_data_2(:,j) = c_data_out_2;   
        end
        [~,~,~,stats] = ttest(shuffled_data_1,shuffled_data_2);
        curr_perm(i,:) = stats.tstat;
    end
    permuted_vals(permi,:,:) = curr_perm;

    
    % saving max for correction for multiple comparisons
    %max_pixel_pvals(permi,:) = [ min(permuted_vals(:)) max(permuted_vals(:)) ];
    
    max_pixel_pvals(permi,:) = [ min(curr_perm(:)) max(curr_perm(:)) ];

    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    % note that here, clusters were obtained by parametrically thresholding
    % the t-maps
    curr_perm(abs(curr_perm)<tinv(1-voxel_pval,n_subs-1))=0;
    clustinfo = bwconncomp(curr_perm);
    max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % notes: cellfun is superfast, and the zero accounts for empty maps
end

%% getting zmaps
zmap = (t_mat-squeeze(mean(permuted_vals,1)))./squeeze(std(permuted_vals));

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
%zmapthresh(zmapthresh>lower_threshold & zmapthresh<upper_threshold)=0; %% change to statistic dist not the z score?
zmapthresh(t_mat>lower_threshold & t_mat<upper_threshold)=0; %% change to statistic dist not the z score?

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



% Check normal distributions

size_of_m = size(permuted_vals);
figure
m = 10;
n = 10;
for i = 1:n*m
 hold on
 data2plot = permuted_vals(:,randi([1 size_of_m(2)]), randi([1 size_of_m(3)]));
 subplot(m,n,i)
 histogram(data2plot)
end

