%% Testing to see what the effect of referencing is on gamma and pre stim alpha


addpath Z:\5_Jessica_Archibald\MRS_ATLAS_2021\Time_frequency_data\Gamma_scaling

%% read ratings
my_subs = {'S12'};
my_extensions = ["_Avg","_Fz_ref","_Mast"];
%% Variables
s_rate = 1000;
time_range = [-500 1000]; % in s
time_vect = linspace(time_range(1),time_range(2),abs(((time_range(1)-time_range(2)))/1000)*s_rate);
n_trials = 20;
n_chans= 31;

%% Looping through the files
for i_subs = 1:length(my_subs)
    S_ID = char(my_subs(i_subs));
    cd(['Z:\30_Oscar_Ortiz\GluEEG\3_Exports\' S_ID])
    my_cd =pwd;
    my_ratings = readtable('Z:\30_Oscar_Ortiz\GluEEG\Ratings.xlsx','Sheet',S_ID);

for i_ext = 1:length(my_extensions)

  

%% Read EEG file and ratings for 1st set
listing = dir(['*_DC Detrend_Normal' char( my_extensions(i_ext)) '*.CSV']);
%% check that there are 4 files
n_files = length(listing);

fprintf(['Processing ' num2str(n_files) ' files for ' S_ID ' with extension ' char(my_extensions(i_ext)) '...'])
fprintf('\n');
for i_n_files = 1:n_files 
[EEGdata1, chan_names1 ] = readEEGfile([my_cd '\' listing(i_n_files).name]);
EEG1 = reshape(EEGdata1,length(time_vect),[],n_chans);
% to clean data
n_trials_recorded_r = height(my_ratings);
% if i_n_files == 3
% %     EEG1(:,end,:) = [];
%     n_trials_recorded_r = 12;
% end
EEG_dims = size(EEG1);
n_trials_recorded = EEG_dims(2);


if n_trials_recorded_r == n_trials_recorded
    fprintf(['Number of EEG segments = number of recorded ratings for ' S_ID ' trial ' listing(i_n_files).name '. Coninuing processing...'])
    fprintf('\n');

    % Generate subject ID structure
    temp1 = my_ratings.Properties.VariableDescriptions{i_n_files};
    temp2 = split(temp1,'_');

    ID_Trial_v = repmat(temp2(1),n_trials_recorded,1);
    ID_Int_v = repmat(temp2(2),n_trials_recorded,1);
    ID_S_v = repmat({S_ID},n_trials_recorded,1);
%     % manual shit
%     ID_Int_v = repmat(temp2(2),20,1);
%     ID_Trial_v = repmat(temp2(1),20,1);
%     ID_S_v = repmat({S_ID},20,1);

    ID_ratings = table2cell(my_ratings(:,i_n_files));
%     if i_n_files == 3
% %       EEG1(:,end,:) = [];
%       ID_ratings = ID_ratings(end-11:end);
%     end
    ID_table1 = table(ID_S_v,ID_Trial_v,ID_Int_v,ID_ratings);
    if i_n_files == 1
        ID_table = ID_table1;
        EEG_out = EEG1;
    else
        ID_table = [ID_table;ID_table1];
        EEG_out = cat(2,EEG_out,EEG1);
    end
    
else
   disp(['Number of EEG segments(' num2str(n_trials_recorded) ')  is NOT number of recorded ratings (20) for ' S_ID ' trial ' listing(i_n_files).name '. Please Check...'])
   fprintf('\n');

end
end



%% Save compiled EEG_variables

save(['Z:\30_Oscar_Ortiz\GluEEG\gamma_sta\' S_ID '_compiled_data' char(my_extensions(i_ext)) '.mat'],"ID_table","EEG_out")

end
end


