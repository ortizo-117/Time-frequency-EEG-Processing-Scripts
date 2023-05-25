function struct_out = wrap_npp_stb_2015(struct_in,names_data,n_permutes,time_window,plotting_input)


struct_out = struct_in;
n_chans = length(struct_in.chanlocs);
n_datasets = length(names_data);
v_time = struct_in.times;
time_s = dsearchn(v_time',time_window(1)); % start time for tf window
time_e = dsearchn(v_time',time_window(2)); % end time for tf window
tfv_time_in  = v_time(time_s:time_e);



for i = 1:n_datasets
    
    x_wait = 0;
    f = waitbar(x_wait,'Performing the permutations...','Name','wrap_npp_stb',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(f,'canceling',0);

    data_in = eval(['struct_in.' char(names_data(i)) '_data;']);
    v_freq = eval(['struct_in.' char(names_data(i)) '_frex;']);
    output_size = zeros(n_chans,length(v_freq),length(tfv_time_in));
    %% initialize storage matrices
    temp.zmap = output_size;
    temp.zmapthresh = output_size;
    temp.zmapthresh_plc = output_size;
    temp.zmapthresh_clc = output_size;
    
    for j = 1:n_chans
        data_chan = squeeze(data_in(:,j,:,:));
        my_title = struct_in.chanlocs(j).labels;
         if getappdata(f,'canceling')
            delete(f)
            msg = 'Error, user terminated the function.';
            error(msg)
         end
         x_wait = j/n_chans;
         waitbar(x_wait,f,['performing the permutations of ' my_title ' in the ' char(names_data(i)) ' range.'] );
    
        [zmap_c, zmapthresh_c, zmapthresh_plc_c, zmapthresh_clc_c, ~] = npp_stb_2015(data_chan,n_permutes,v_time,v_freq,time_window,plotting_input,my_title);
        temp.zmap(j,:,:) = zmap_c;
        temp.zmapthresh(j,:,:) = zmapthresh_c;
        temp.zmapthresh_plc(j,:,:) = zmapthresh_plc_c;
        temp.zmapthresh_clc(j,:,:) = zmapthresh_clc_c;
    end
    eval(['struct_out.' char(names_data(i)) '_zmaps = temp;']);
    delete(f)
end


end
