% A script that uses EEGLAB's computational engine to process .set files. Modified by Oscar
% Processing includes: 
% 1. Importing Events from Digi
% 2. Importing electrode Locations from file
% 3. Filtering highpass (0.5 Hz), lowpass (40 Hz)
% 4. Epoching based on events
% 5. Rereferencing to average
% 6. Displaying Topoplot of the epoch average for each electrode
% 7. Saves (overwrites) the files in the folder where they are located. 
% Inputs: pathname - the directory location of the folder holding the .set
% and .fdt files you want to process. Eg. 'Z:\19_Carson_Berry\EEG':
% Outputs: Saved files, topoplots. 
% Files to be check manually for noise after processing and before ICA. 



function preproces_EEG_folder(pathname,saving_folder,extension_to_look_for)
 avg_count=0; %init avg_count that counts the instances where both mastoids are not available. 
 ref_boolean = 0;
 file_struct_list = dir([pathname filesep() extension_to_look_for]);  %% get list of .set files in the pathname specified

 filename_cell_list = {file_struct_list.name};  %% extract the filenames into a cellarray

 filename_list=deblank(char(filename_cell_list));

 endout=regexp(pathname,filesep,'split');
 saving_name = [char(endout(end-1)) '_' char(endout(end))];

 length_filename=size(filename_list);
 [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;  
 
 for k = 1:length_filename(1)
          %subject = file_struct_list.name  %% this iterates over the elements of the cell array, one-by-one, setting the `filename` variable like a loop variable
          % save(debugging_variables);
          
          filename=deblank(filename_list(k, :));
%           filename = erase(filename,'.poly5');
%           filename = [filename '.set'];
          
          [~,filename_text,~]= fileparts(strcat(pathname,'\',filename));
          
          EEG.etc.eeglabvers = '14.1.1'; % this tracks which version of EEGLAB is being used, you may ignore it
          EEG = pop_loadset('filename',filename,'filepath',pathname); %loads the specified file into eeglab
          
          [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', 'raw');
          EEG = eeg_checkset( EEG );
          %eeg.nbchan-1 is the second last eeg channel, digi.
          
          EEG = pop_chanevent(EEG, EEG.nbchan-1,'oper','X<255','edge','leading','edgelen',1,'duration','on','delchan','off','delevent', 'on','nbtype',1);
          EEG = eeg_checkset( EEG );
          Event_chan=strcat('chan',int2str(EEG.nbchan-1)); 
          EEG = pop_select( EEG,'nochannel',{'Digi' 'Saw'});  %remove digi and saw channels
          EEG.setname = strcat(filename_text, '_no_digi');          %sets name to filename with '_no_digi' appended to end
          EEG = eeg_checkset( EEG );
          
          try
              fprintf('Looking up (%i) channels in %s', EEG.nbchan, filename);
              mroot=matlabroot; %get matlab root
              EEG=pop_chanedit(EEG, 'lookup',[mroot '\\toolbox\\eeglab14_1_1b\\plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp']); % get cap locations from matlabroot folder
              EEG = eeg_checkset( EEG );
              
          catch
               warning('EEG channel names not provided. Looking up from archive electrode location files.');
              %this line will work on any kramer computer networked to the z: drive.
              if(EEG.nbchan==29)
                  EEG = pop_chanedit(EEG, 'load',{'Z:\19_Carson_Berry\EEG\MATLAB\trunk\src\cap locations\29 channel, missing f3, m1, m2.ced' 'filetype' 'autodetect'}); %if this program is being run on a MAC- this location file path will need to use \\ instead of \
                  EEG = eeg_checkset( EEG );
              else
                  if(EEG.nbchan==31)
                      EEG = pop_chanedit(EEG, 'load',{'Z:\19_Carson_Berry\EEG\MATLAB\trunk\src\cap locations\31_channel_locs_missing_F3.ced' 'filetype' 'autodetect'}); %if this program is being run on a MAC- this location file path will need to use \\ instead of \
                      EEG = eeg_checkset( EEG );
                  else
                      if(EEG.nbchan==32)
                          EEG = pop_chanedit(EEG, 'load',{'Z:\19_Carson_Berry\EEG\MATLAB\trunk\src\cap locations\32_channel_locs_noice.ced' 'filetype' 'autodetect'}); %if this program is being run on a MAC- this location file path will need to use \\ instead of \
                          EEG = eeg_checkset( EEG );
                      else
                      end
                  end
              end
          end

          %EEG = pop_chanedit(EEG, 'lookup','C:\Users\Kip\Documents\eeglab2021.1\plugins\dipfit\standard_BESA\standard-10-5-cap385.elp');
          EEG = pop_resample( EEG, 500);
          EEG = eeg_checkset( EEG );
          EEG = pop_eegfiltnew(EEG, 'locutoff',0.1,'hicutoff',120);          %EEG.setname=strcat(filename_text,'_highpass');
          EEG = eeg_checkset( EEG );
          EEG = pop_eegfiltnew(EEG, 'locutoff',59,'hicutoff',61,'revfilt',1);
          EEG = eeg_checkset( EEG );
          
          num_chan=EEG.nbchan; % get number of channels
          % To plot the N1 latencies we will rereference to FZ
          chan_labels = {EEG.chanlocs.labels}.'; %make a simple cell array out of the cell array within the EEG.chanlocs struct
                
                          %To plot the N2P2 latencies we will reference to m1/m2
             M1_index= find(strcmpi(chan_labels, 'M1'));
             M2_index= find(strcmpi(chan_labels, 'M2'));                                

             if((isempty(M1_index)==0)&&(isempty(M2_index)==0))
                 %referencing to M1 and M2
                 EEG = pop_reref( EEG, [M1_index M2_index] ,'keepref','on'); %rereference to M1/M2
                 [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
                 EEG = eeg_checkset( EEG );
                 EEG.setname=strcat(filename_text,'_Ear_ref');
                 ref_boolean = 1 ; % if there is mastoids, make the ref boolean 1

                 
                 
             else
                 %referencing to average because M1/M2 are not
                 %recorded.
                 
                 EEG = pop_reref( EEG, []); %rereference to average
                 [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
                 EEG = eeg_checkset( EEG );
                 EEG.setname=strcat(filename_text,'_Avg_ref');
                 ref_boolean = 2; % if we use avg reference make the boolean = 2;
           
             end
             
                 
                 
          EEG = eeg_checkset( EEG ); 
          EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on','pca',26);
          EEG = eeg_checkset( EEG ); 

           EEG = pop_epoch( EEG, {  Event_chan  }, [-1  1], 'newname', strcat(filename_text,'_epoched'), 'epochinfo', 'yes');
           EEG = eeg_checkset( EEG );
           EEG = pop_rmbase( EEG, [-500     0]);
           EEG.setname=strcat(filename_text,'_reref');
          
          % saving
          if ref_boolean == 1 % if we used the mastoids
              ref_ext = '_Ear_ref';
          elseif ref_boolean == 2 % if we used the average reference
              ref_ext = '_Avg_ref';
          else
              ref_ext = '_wtf_happened';
          end
          
          saving_name = [filename_text ref_ext '_filt_ICA_epoch'];
          EEG = pop_saveset( EEG, 'filename',saving_name,'filepath', saving_folder);
          EEG = eeg_checkset( EEG );
          [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
          eeglab redraw;
          fprintf('\n\n\n %i percent done processing folder \n\n\n',k/length_filename(1)*100);         
 end
 
 fprintf('All done Processing!');
 
end

