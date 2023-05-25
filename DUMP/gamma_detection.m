

function files_not_processed  = gamma_detection(pathname)
%Author: Cassie Choles
%Date: Oct 29/21 
%Description: This function is to plot gamma activity
%The input of this function is the directory where the .set files are,
%that have already been decomposed by ICA.
fnp_i = 1;
problem_v = [];
file_struct_list = dir([pathname filesep() '*.set']);  %% get list of .set files in the pathname specified

 filename_cell_list = {file_struct_list.name};  %% extract the filenames into a cellarray

 filename_list=deblank(char(filename_cell_list));

 endout=regexp(pathname,filesep,'split');
 saving_name = [char(endout(end-1)) '_' char(endout(end))];

 length_filename=size(filename_list);
 
 [ALLEEG EEG CURRENTSET ALLCOM] = eeglab  
 
 
  for k = 1:length_filename(1)
     %% load file
     filename=deblank(filename_list(k, :));   
     [~,filename_text,~]= fileparts(strcat(pathname,'\',filename));   
%      EEG.etc.eeglabvers = '14.1.1'; % this tracks which version of EEGLAB is being used, you may ignore it
     EEG = pop_loadset('filename',filename,'filepath',pathname); %loads the specified file into eeglab
     
 
     
     %%  extracting gamma for each channel
     try
     num_chan=EEG.nbchan; % get number of channels
     gamma_mat = zeros(num_chan,19,400);
     lf_mat = zeros(num_chan,8,400);
     for i = 1:num_chan
        try
        fh1 = figure;
        [ersp_gamma, ~, ~, times_gamma, frequencies_gamma, erspboot, ~, ~] = pop_newtimef( EEG, 1, i , [-500  999], [0] , 'topovec', i, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', char(EEG.chanlocs(i).labels), 'baseline',[-499 -100], 'freqs', [30 100], 'freqscale', 'linear', 'plotphase', 'off', 'ntimesout', 400, 'padratio', 1);
        gamma_mat(i,:,:) = ersp_gamma;
        pause(1)
        close(fh1)
        fh2 = figure;
        [ersp_lf, ~, ~, times_lf, frequencies_lf, erspboot, ~, ~] = pop_newtimef( EEG, 1, i , [-500  999], [0] , 'topovec', i, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', char(EEG.chanlocs(i).labels), 'baseline',[-499 -100], 'freqs', [0.5 30], 'freqscale', 'linear', 'plotphase', 'off', 'ntimesout', 400, 'padratio', 1);
        lf_mat(i,:,:) = ersp_lf;
        pause(1)
        close(fh2)
        catch
        disp(['Problem with channel '  char(EEG.chanlocs(i).labels) '. Potentially all zeros. Pressing on...'])
        problem_v = [problem_v , i];
        end       
     end
     %% put it into EEG structure 
     
     EEG.gamma_data = gamma_mat;
     EEG.gamma_frex = frequencies_gamma;
     EEG.gamma_time = times_gamma;
     EEG.lf_data = lf_mat;
     EEG.lf_frex = frequencies_lf;
     EEG.lf_time = times_lf;  
     EEG.chans_with_problems = problem_v;
     
     
     %% Save file
     saving_name = [filename_text '_tf_extracted'];
     EEG = pop_saveset( EEG, 'filename',saving_name,'filepath', pathname);
     EEG = eeg_checkset( EEG );
     [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
     eeglab redraw;
     fprintf('\n\n\n %i percent processing folder \n\n\n',k/length_filename(1)*100); 
     catch
         fprintf([filename 'not processed. Please check file']);
         if fnp_i == 1
             files_not_processed = string(filename);
             fnp_i = fnp_i+1;
         else
             files_not_processed =vertcat(files_not_processed,string(filename));
             fnp_i = fnp_i+1;
         end
     end
     
         
  end
  
 

end

