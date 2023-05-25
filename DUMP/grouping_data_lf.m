function struct_out = grouping_data_lf(pathname)

% Groups all gamma data in the folder in one condition


 file_struct_list = dir([pathname filesep() '*low_freq.set']);  %% get list of .set files in the pathname specified
 filename_cell_list = {file_struct_list.name};  %% extract the filenames into a cellarray
 filename_list=deblank(char(filename_cell_list));
 %endout=regexp(pathname,filesep,'split');
 length_filename=size(filename_list);
 
 
 [ALLEEG EEG CURRENTSET ALLCOM] = eeglab ; 
 L_handed_counter = 0;
 
 for k = 1:length_filename(1)
     filename=deblank(filename_list(k, :));   
     [~,filename_text,~]= fileparts(strcat(pathname,'\',filename));   
     EEG = pop_loadset('filename',filename,'filepath',pathname); %loads the specified file into eeglab
     
     % on the first iteration, create output mat with appropriate dimensions
     if k == 1
         output_dims = size(EEG.C3_lf);
         Cc_output = zeros(length_filename(1),output_dims(1),output_dims(2));
         Ci_output = zeros(length_filename(1),output_dims(1),output_dims(2));
         
     else
     end
     
     % add data to output matrix
     if strcmp(filename(1:15),'Participant_C15')
         Cc_output(k,:,:) =  EEG.C4_lf;
         Ci_output(k,:,:) =  EEG.C3_lf;  
         L_handed_counter = L_handed_counter +1;
     else
         Cc_output(k,:,:) =  EEG.C3_lf;
         Ci_output(k,:,:) =  EEG.C4_lf;
     end  
 end
 
 struct_out.Cc_lf = Cc_output;
 struct_out.Cc_mean = squeeze(mean(Cc_output));
 struct_out.Ci_lf = Ci_output;
 struct_out.Ci_mean = squeeze(mean(Ci_output));
 struct_out.times = EEG.tf_time;
 struct_out.frex = EEG.tf_frex;
 struct_out.filenamelist = filename_list;
 struct_out.L_handed_counter = L_handed_counter;

      
 end
 
