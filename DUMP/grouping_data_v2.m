function struct_out = grouping_data_v2(pathname)

% Groups all gamma data in the folder in one condition


 file_struct_list = dir([pathname filesep() '*.set']);  %% get list of .set files in the pathname specified
 filename_cell_list = {file_struct_list.name};  %% extract the filenames into a cellarray
%  B = cellfun(@(x) strsplit(x,'_'),filename_cell_list,'UniformOutput',false);
%  mask = zeros(length(B),1);
%  for i = 1 : length(B)
%      if string(B{1,i}(3))== area_to_subset && string(B{1,i}(1))== timepoint
%          mask(i) = 1;
%      end
%  end
%  
%  new_filename_cell_list = filename_cell_list(logical(mask));
 filename_list=deblank(char(filename_cell_list));
 
 length_filename=size(filename_list);
 
 
 [ALLEEG EEG CURRENTSET ALLCOM] = eeglab ; 
 L_handed_counter = 0;
 
 for k = 1:length_filename(1)
     filename=deblank(filename_list(k, :));   
     [~,filename_text,~]= fileparts(strcat(pathname,'\',filename));   
     EEG = pop_loadset('filename',filename,'filepath',pathname); %loads the specified file into eeglab
     
     % on the first iteration, create output mat with appropriate dimensions
     if k == 1
         output_dims = size(EEG.C3_gamma);
         Cc_output = zeros(length_filename(1),output_dims(1),output_dims(2));
         Ci_output = zeros(length_filename(1),output_dims(1),output_dims(2));
         
     else
     end
     
     % add data to output matrix
     if strcmp(filename(1:15),'Participant_C15')
         Cc_output(k,:,:) =  EEG.C4_gamma;
         Ci_output(k,:,:) =  EEG.C3_gamma;  
         L_handed_counter = L_handed_counter +1;
     else
         Cc_output(k,:,:) =  EEG.C3_gamma;
         Ci_output(k,:,:) =  EEG.C4_gamma;
     end  
 end
 
 struct_out.Cc_gamma = Cc_output;
 struct_out.Cc_mean = squeeze(mean(Cc_output));
 struct_out.Ci_gamma = Ci_output;
 struct_out.Ci_mean = squeeze(mean(Ci_output));
 struct_out.times = EEG.tf_time;
 struct_out.frex = EEG.tf_frex;
 struct_out.filenamelist = filename_list;
 struct_out.L_handed_counter = L_handed_counter;

      
 end