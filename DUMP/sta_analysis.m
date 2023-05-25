

% load eeglab
addpath  C:\Users\User\OneDrive\Documents\MATLAB\eeglab2021.0

subs2process = {'S23','S24','S25','S26','S27','S28','S29'};
my_filepath = 'Z:\30_Oscar_Ortiz\GluEEG\gamma_sta';
save_name = '_saved_set.set';

for i = 1:length(subs2process)
    load([my_filepath filesep() char(subs2process(i)) '_tf_data__Mast.mat'])
    EEG = pop_saveset( EEG, 'filename',[char(subs2process(i)), save_name],'filepath',my_filepath);
    
end

