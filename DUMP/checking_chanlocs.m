
function  checking_chanlocs(pathname)

file_struct_list = dir([pathname filesep() '*.set']);  %% get list of .set files in the pathname specified

 filename_cell_list = {file_struct_list.name};  %% extract the filenames into a cellarray

 filename_list=deblank(char(filename_cell_list));

 endout=regexp(pathname,filesep,'split');
 saving_name = [char(endout(end-1)) '_' char(endout(end))];

 length_filename=size(filename_list);
 
 [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;  
 
 for k = 1:length_filename(1)
     
     filename=deblank(filename_list(k, :));
%    filename = erase(filename,'.poly5');
%    filename = [filename '.set'];
          
      [~,filename_text,~]= fileparts(strcat(pathname,'\',filename));
          
      EEG.etc.eeglabvers = '14.1.1'; % this tracks which version of EEGLAB is being used, you may ignore it
      EEG = pop_loadset('filename',filename,'filepath',pathname); %loads the specified file into eeglab
      
      c_x = EEG.chanlocs(1).X;
      c_z = EEG.chanlocs(1).Z;
      c_y = EEG.chanlocs(1).Y;
      
      u_x = EEG.chanlocs(7).X;
      u_z = EEG.chanlocs(7).Z;
      u_y = EEG.chanlocs(7).Y;
      
      if c_x == u_x && c_z == u_z && c_y == u_y
          EEG=pop_chanedit(EEG, 'load',[],'load',[],'load',{'Z:\\19_Carson_Berry\\EEG\\MATLAB\\trunk\\src\\cap locations\\31_channel_locs_missing_F3.ced' 'filetype' 'autodetect'});
            
      end
      % saving
      saving_name = [filename_text '_chanlocs_fixed'];
      EEG = pop_saveset( EEG, 'filename',saving_name filepath', pathname);
      EEG = eeg_checkset( EEG );
      [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
      eeglab redraw;
      fprintf('\n\n\n %i percent done processing folder \n\n\n',k/length_filename(1)*100);        
 
 
 end
 
 


end

