

function gamma_detection_no_baseline(pathname,saving_folder)
%Author: Cassie Choles
%Date: Oct 29/21 
%Description: This function is to plot gamma activity
%The input of this function is the directory where the .set files are,
%that have already been decomposed by ICA.

file_struct_list = dir([pathname filesep() '*Ear_ref_filt_ICA_epoch_pruned.set']);  %% get list of .set files in the pathname specified

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
     
     %% Epoching  
%      EEG = pop_epoch( EEG, {  'chan33'  }, [-.4  1], 'newname', 'Participant_C07_condition_2_Ear_ref epochs', 'epochinfo', 'yes');
%      EEG = eeg_checkset( EEG );
%      EEG = pop_rmbase( EEG, [-400 0] ,[]);
%      EEG = eeg_checkset( EEG );

     
     %%  finding the index of C3 & C4 
     num_chan=EEG.nbchan; % get number of channels
     chan_labels = {EEG.chanlocs.labels}.'; %make a simple cell array out of the cell array within the EEG.chanlocs struct
         
     C3_index= find(strcmpi(chan_labels, 'C3'));
     C4_index= find(strcmpi(chan_labels, 'C4'));                                

     
     %% extract gamma 
     fh1 = figure;
     [ersp, ~, ~, times, frequencies, ~, ~, ~] = pop_newtimef( EEG, 1, C3_index, [-400  999], [0] , 'topovec', C3_index, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'C3', 'baseline',[NaN], 'freqs', [30 100], 'freqscale', 'log', 'plotphase', 'off', 'ntimesout', 400, 'padratio', 1);
     C3_gamma = ersp;
     pause(1)
     
     %[ersp, ~, ~, times, frequencies, erspboot, ~, ~] = pop_newtimef( EEG, 1, 14, [-200  999], [0] , 'topovec', 14, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'C3', 'baseline',[NaN], 'alpha',0.05, 'freqs', [30 100], 'plotphase', 'off', 'padratio', 1);

     
     fh2 = figure;
     [ersp, ~, ~, times, frequencies, ~, ~, ~] = pop_newtimef( EEG, 1, C4_index, [-400  999], [0] , 'topovec', C4_index, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'C4', 'baseline',[NaN], 'freqs', [30 100], 'freqscale', 'log', 'plotphase', 'off', 'ntimesout', 400, 'padratio', 1);
     C4_gamma = ersp;
     tf_frex = frequencies;
     tf_time = times;
     pause(1)
     close(fh2)
     close(fh1)
     
     %% put it into EEG structure 
     
     EEG.C3_gamma = C3_gamma;
     EEG.C4_gamma = C4_gamma;
     EEG.tf_frex = tf_frex;
     EEG.tf_time = tf_time;
     
     
     %% Save file
     saving_name = [filename_text '_gamma'];
     EEG = pop_saveset( EEG, 'filename',saving_name,'filepath', saving_folder);
     EEG = eeg_checkset( EEG );
     [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
     eeglab redraw;
     fprintf('\n\n\n %i percent processing folder \n\n\n',k/length_filename(1)*100); 
     
  end
  
 

end
