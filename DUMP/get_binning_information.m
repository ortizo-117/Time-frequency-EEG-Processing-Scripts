% Author: Oscar Ortiz
% Date: 2021/05/19
% Purpose: Read in raiting data from each participant, scales the raitings
%   from 1 to 100, then bins the individual trials into four categories based
%   on the scaled raiting (BIN 1 = 0 to 25, BIN 2 = 25 to 50, BIN 3 = 50 to
%   75, BIN 4 = 75 to 100).Each bin includes the left bin edge, except for 
%   the last bin, which includes both bin edges.
% Inputs: string of the file where the data is stored in a 20 by 8 fashion
%   (ie. 1st colum is the trial count and second one is the actual raiting
%   and so on).
% Outputs = matrix containing the binning index per trial.

function binning_index = get_binning_information(file_name_in)

% Read in xlxs as table
t = readtable(file_name_in);
datain = t{:,:};

% Get the raitings from the data only (every second column)
ratings = datain(:,2:2:end);
raitings_long = ratings(:);

% Scale raitings from 0 to 100
low_end = min(raitings_long); high_end =max(raitings_long); 
scaled_raitings_long = 0 + (((100 - 0) * (raitings_long - low_end))/(high_end - low_end));

% Bin the data
binning_index = discretize(scaled_raitings_long,[0 25 50 75 100]); % here are the bin limits

% % Reshaping output
% binned_data = reshape(binning_index,[20 4]);

end