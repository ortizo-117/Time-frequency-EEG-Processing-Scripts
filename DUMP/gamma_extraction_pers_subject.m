%% Gamma/N2P2 extraction based on rating 
addpath C:\Users\ortizo\Documents\MATLAB\eeglab2021.1
eeglab
cd Z:\30_Oscar_Ortiz\GluEEG\gamma_sta
s_list = ["S12"];
ref_list = ["_Mast"];
% subs missing S03 (justified), S12,
for s_i = 1: length(s_list)
    c_sub = char(s_list(s_i));
    for r_i = 1:length(ref_list)
        c_ref = char(ref_list(r_i));
       
load([c_sub '_compiled_data' c_ref '.mat'])

% get a subset of trials with only 3 or 4
T = ID_table;
my_mask = cell2mat(T.ID_ratings) == 3 | cell2mat(T.ID_ratings) == 4;
trials2get = find(my_mask==1);

% if there is more than 20, only get 20 random trials
n_stim_sub = sum(my_mask);



if n_stim_sub >= 20
    trials2get = trials2get(randperm(20));
elseif n_stim_sub < 20 && n_stim_sub >=11
    trials2get = trials2get;
elseif n_stim_sub < 11
    trials2get = [];
end

EEG_subset = EEG_out(:,trials2get,:);
if ~isempty(EEG_subset)

%% Create an EEG structure to trick EEG lab into processing data
load('chan_info.mat')

% creating necessary variables
s_rate = 1000;
time_range = [-500 1000]; % in s
time_vect = linspace(time_range(1),time_range(2),abs(((time_range(1)-time_range(2)))/1000)*s_rate);
n_trials = size(EEG_subset,2);
n_chans= size(EEG_subset,3);
pnts = size(EEG_subset,1);

% reorganize EEG data
EEG_reshaped = permute(EEG_subset,[3 1 2]);

%creating EEGlab Structure
EEG.data= EEG_reshaped;
EEG.srate = s_rate;
EEG.times = time_vect;
EEG.chanlocs = chanlocs;
EEG.chaninfo = chaninfo;
EEG.xmin = min(EEG.times)/s_rate;
EEG.xmax = max(EEG.times)/s_rate;
EEG.pnts = length(time_vect);
EEG.nbchan = size(EEG_reshaped,1);


%% Gamma detection 


epoch_length = [-500 1000];
freqs_in = [30 100];
baseline_in = [-450 -50];
n_times_out = 1000;
pads = 4 ;% can be 2^n 
win_size = 150;
NF = [];
c_lims = [-4 4];

EEG_out = gamma_detection_sta(EEG,epoch_length,freqs_in,baseline_in,n_times_out,pads,win_size,NF,c_lims);
EEG = EEG_out;

save([c_sub '_tf_data_' c_ref],'EEG');
else 
    fprintf([c_sub c_ref ' does not have any trials. Not processing...'])

end

    end
end



