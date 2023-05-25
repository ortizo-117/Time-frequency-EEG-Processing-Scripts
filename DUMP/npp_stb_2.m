function [zmap, zmapthresh, zmapthresh_plc, zmapthresh_clc, tfv_time, v_freq] = npp_stb_2(input_mat,n_permutes,v_time,v_freq,time_window,baseline_range,plotting_input,my_title)
%% Description:
 % Takes a time frequency matrix (dimensions = n_subs, frequencies, times)
 % and peforms the non perametric permutation comparing the signal to the
 % baseline. Inputs include input mat, number of permutations (n_permutes),
 % vectors for time and frequency(v_time and v_freq), limits in the time 
 % axis where you want to perform the analysis (time_window), and the
 % plotting_input variable, in which if is set to 1, it will plot the 
 % responses. The output consist of the z_map, the logical mask for the 
 % zmap uncorrected for multple comparisons (zmapthresh), the logical mask
 % with a pixel level correction for multiple comparisons (zmapthresh_plc),
 % and the logical threshold map with cluster correction for mutliple
 % comparisons (zmapthesh_clc). The last output is the new time vector
 % based on time_window. 



%% getting parameters for permutations
time_s = dsearchn(v_time',baseline_range(1)); % start time for tf window
time_0 = dsearchn(v_time',baseline_range(2));
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
tfv_time  = v_time(time_s:time_e);
nTimepoints = numel(tfv_time);
avg_mat =squeeze(mean(input_mat));
temp = size(input_mat);
n_subs = temp(1);

% statistical significant values for voxel and cluster analysis
voxel_pval   = 0.005;
cluster_pval = 0.005;
mcc_voxel_pval = 0.05; % mcc = multiple comparisons correction
mcc_cluster_pval = 0.0001;
%% Running actual T test on data
% need to run it per frequency bin

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

%% running the premutations

% initialize null hypothesis matrices
permuted_maxvals = zeros(n_permutes,2,length(v_freq));
max_clust_info   = zeros(n_permutes,1);
max_pixel_pvals = zeros(n_permutes,2);
%% permutation loop 
curr_perm = zeros(mat_size); % initializing a variable to temporary store the results of current permutation before reshaping back into tf spectograph
permuted_vals = zeros([n_permutes, mat_size]);

% waitbar
x_wait = 0;
f = waitbar(x_wait,'Performing the permutations...','Name','npp_stb',...
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
    for i = 1:mat_size(1) % per frequency bin 
        prestim_estimate = squeeze(input_mat(:,i,time_s:time_0));
        poststim_points = squeeze(input_mat(:,i,time_s:time_e));
        pre_stim_estimate_r = zeros([n_subs,1]);
        temp1 = size(prestim_estimate);
        for m = 1:n_subs
            r_index = randi([1,temp1(2)]);
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
    max_pixel_pvals(permi,:) = [ min(curr_perm(:)) max(curr_perm(:)) ];

    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    % note that here, clusters were obtained by parametrically thresholding
    % the t-maps
    curr_perm(abs(curr_perm)<tinv(1-voxel_pval,n_subs-1))=0;
    clustinfo = bwconncomp(curr_perm);
    max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % notes: cellfun is superfast, and the zero accounts for empty maps
end
delete(f)

%% getting zmaps
zmap = (t_mat-squeeze(mean(permuted_vals,1)))./squeeze(std(permuted_vals));
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
my_clims = [-9 9];
if plotting_input == 1
    figure
    my_lxlim = time_window(1)+10;
    my_hxlim = time_window(2);
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