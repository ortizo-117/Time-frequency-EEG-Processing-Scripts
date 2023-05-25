function=
%Author: Cassie Choles
%Date: Oct 29/21 
%Description: This function is to plot gamma activity
%The input of this function is the directory where the .set files are,
%that have already been decomposed by ICA.
%addpath C:\Users\Kip\Documents\eeglab2021.1 % add eeglab path
%% Indicate path with all folders 

cd  Z:\30_Oscar_Ortiz\LEPs_CHEPs\scripts_functions
addpath Z:\30_Oscar_Ortiz\LEPs_CHEPs\Raw
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\Raw\processed_data'; % folder 
addpath C:\Users\Kip\Documents\eeglab2021.1


file_struct_list = dir([pathname filesep() 'epoch_chanlocs_fixed_pruned.set']);  %% get list of .set files in the pathname specified

 filename_cell_list = {file_struct_list.name};  %% extract the filenames into a cellarray

 filename_list=deblank(char(filename_cell_list));

 endout=regexp(pathname,filesep,'split');
 saving_name = [char(endout(end-1)) '_' char(endout(end))];

 length_filename=size(filename_list);
 
 [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;  
  for k = 1:length_filename(1)
     filename=deblank(filename_list(k, :));
     %for electrodes C3 & C4 electrodes:
     EEG.etc.eeglabvers = '14.1.1'; % this tracks which version of EEGLAB is being used, you may ignore it
     EEG.etc.eeglabvers = '2021.1'; % this tracks which version of EEGLAB is being used, you may ignore it
     EEG.setname='raw';
     EEG = eeg_checkset( EEG );
     EEG = eeg_checkset( EEG );
     EEG.etc.eeglabvers = '14.1.1'; % this tracks which version of EEGLAB is being used, you may ignore it
     EEG.etc.eeglabvers = '2021.1'; % this tracks which version of EEGLAB is being used, you may ignore it
     EEG = eeg_checkset( EEG );
     EEG = eeg_checkset( EEG );
     figure; output = pop_newtimef( EEG, 1, 15, [-200  999], [0] , 'topovec', 15, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'C3', 'baseline',[0], 'freqs', [30 100], 'freqscale', 'log', 'plotphase', 'off', 'ntimesout', 400, 'padratio', 1);
     figure; output = pop_newtimef( EEG, 1, 17, [-200  999], [0] , 'topovec', 17, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'C4', 'baseline',[0], 'freqs', [30 100], 'freqscale', 'log', 'plotphase', 'off', 'ntimesout', 400, 'padratio', 1);
     saving_name = [filename_text '_gamma'];
     EEG = pop_saveset( EEG, 'filename',saving_name,'filepath', pathname);
     EEG = eeg_checkset( EEG );
     [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
     eeglab redraw;
     fprintf('\n\n\n %i percent processing folder \n\n\n',k/length_filename(1)*100); 
     
  end
  
 

 