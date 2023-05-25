function struct_out = wrap_npp_c1tc2(struct_in1,struct_in2,names_data,n_permutes,time_window,plotting_input)



n_chans_1 = length(struct_in1.chanlocs);
n_chans_2 = length(struct_in1.chanlocs);
if n_chans_1 ~= n_chans_2
    msg = 'Error, trying to compare two datasets with different number of channels';
    error(msg)
end

n_datasets = length(names_data);
v_time = struct_in1.times;
time_s = dsearchn(v_time',time_window(1)); % start time for tf window
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
tfv_time_in  = v_time(time_s:time_e);



for i = 1:n_datasets
    
    x_wait = 0;
    f = waitbar(x_wait,'Performing the permutations...','Name','wrap_npp_c1tc2',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(f,'canceling',0);

    data_in_1 = eval(['struct_in1.' char(names_data(i)) '_data;']);
    data_in_2 = eval(['struct_in2.' char(names_data(i)) '_data;']);

    v_freq = eval(['struct_in1.' char(names_data(i)) '_frex;']);
    output_size = zeros(n_chans_1,length(v_freq),length(tfv_time_in));
    
    %% initialize storage matrices
    temp.zmap = output_size;
    temp.zmapthresh = output_size;
    temp.zmapthresh_plc = output_size;
    temp.zmapthresh_clc = output_size;
    
    for j = 1:n_chans_1
        data_chan_1 = squeeze(data_in_1(:,j,:,:));
        data_chan_2 = squeeze(data_in_2(:,j,:,:));
        my_title = struct_in2.chanlocs(j).labels;
         if getappdata(f,'canceling')
            delete(f)
            msg = 'Error, user terminated the function.';
            error(msg)
         end
         x_wait = j/n_chans_1;
         waitbar(x_wait,f,['performing the permutations of ' my_title ' in the ' char(names_data(i)) ' range.'] );
    
        [zmap_c, zmapthresh_c, zmapthresh_plc_c, zmapthresh_clc_c,...
            tfv_time_c, ~] = npp_c1tc2(data_chan_1,data_chan_2,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);
        
        temp.zmap(j,:,:) = zmap_c;
        temp.zmapthresh(j,:,:) = zmapthresh_c;
        temp.zmapthresh_plc(j,:,:) = zmapthresh_plc_c;
        temp.zmapthresh_clc(j,:,:) = zmapthresh_clc_c;
        temp.new_time = tfv_time_c;
    end
    
    
    eval(['struct_out.' char(names_data(i)) '_zmaps = temp;']);
    delete(f)
end


end
