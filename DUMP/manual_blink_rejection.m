



function  manual_blink_rejection( pathname )
%Author: Cassie Choles
%Date: Oct 14/21 
%Description: This function is to manually clean EEG data by visual 
%inspection of ICA components and checking for trials.
%The input of this function is the directory where the .set files are,
%that have already been decomposed by ICA.
%User has to inpsect the components visually to determine which ones 
%to keep.

%% blink and artifact detection step
addpath C:\Users\Kip\Documents\eeglab2021.1
pathname = 'Z:\30_Oscar_Ortiz\LEPs_CHEPs\preprocessed_data';
%manual_blink_rejection(pathname);
file_struct_list = dir([pathname filesep() '*ICA_epoch.set']);  %% get list of .set files in the pathname specified

 filename_cell_list = {file_struct_list.name};  %% extract the filenames into a cellarray

 filename_list=deblank(char(filename_cell_list));

 endout=regexp(pathname,filesep,'split');
 saving_name = [char(endout(end-1)) '_' char(endout(end))];

 length_filename=size(filename_list);
 
 [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;  
 
 for k = 1:length_filename(1)
     filename=deblank(filename_list(k, :));   
     [~,filename_text,~]= fileparts(strcat(pathname,'\',filename));   
      EEG.etc.eeglabvers = '14.1.1'; % this tracks which version of EEGLAB is being used, you may ignore it
      EEG = pop_loadset('filename',filename,'filepath',pathname); %loads the specified file into eeglab
      %Inspsection of components by map
      EEG.etc.eeglabvers = '2021.1'; % this tracks which version of EEGLAB is being used, you may ignore it
     
     j=0;
     while j == 0 % in case you mess up a component selection you can try again
     % Inspect components reflecting blinks and noise visually
     pop_selectcomps(EEG, [1:7] );
     pause 
     prompt = 'Are you happy with your selection of components? (input 1 to continue or 0 to do again)';
     j = input(prompt);
     end
    
%      pop_selectcomps(EEG,[1:7] );
%      pause
%       prompt = 'Blink Component selection ok? (1 = yes, continue. 0 = no, do again)';
%       x = input(prompt);
%       if x == 0
%          pop_selectcomps(EEG, [1:7] );
%          x = input(prompt);
%       else
%       end  
      % pause
       EEG = pop_subcomp( EEG,[], 0);
       EEG = eeg_checkset( EEG );
       
       m=0;
       while m == 0;
         pop_eegplot( EEG, 0, 1, 1);
         pause 
         prompt = 'Are you happy with your deletion of segments? (input 1 to continue or 0 to do again)';
         m = input(prompt);
       end
       
       
        saving_name = [filename_text '_pruned'];
        EEG = pop_saveset( EEG, 'filename',saving_name,'filepath', pathname);
        EEG = eeg_checkset( EEG );
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        eeglab redraw;
        fprintf('\n\n\n %i percent processing folder \n\n\n',k/length_filename(1)*100);      
 end
end

