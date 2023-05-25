

function filter_files(pathname,saving_folder,extension_to_look_for,saving_extension,filter_ranges)
 
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
          EEG = pop_eegfiltnew(EEG, 'locutoff',filter_ranges(1),'hicutoff',filter_ranges(2));          %EEG.setname=strcat(filename_text,'_highpass');
          saving_name = [filename_text  saving_extension];
          EEG = pop_saveset( EEG, 'filename',saving_name,'filepath', saving_folder);
          EEG = eeg_checkset( EEG );
          [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
          eeglab redraw;
          fprintf('\n\n\n %i percent done processing folder \n\n\n',k/length_filename(1)*100);         
 end
 
 fprintf('All done Processing!');
 
end

