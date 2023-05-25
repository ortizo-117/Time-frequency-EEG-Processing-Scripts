% Author: Oscar Ortiz
% Date: 2021/05/21
% Purpose: Bins the individual trials into four categories based
%   on the scaled raiting (BIN 1 = 0 to 25, BIN 2 = 25 to 50, BIN 3 = 50 to
%   75, BIN 4 = 75 to 100).Each bin includes the left bin edge, except for 
%   the last bin, which includes both bin edges.
% Inputs: vector with raitings
% Outputs: vector with binned raitings

function output = bin_raitings_single_chan(raitings)

% Scale raitings from 0 to 100
low_end = min(raitings); high_end =max(raitings); 
scaled_raitings = 0 + (((100 - 0) * (raitings - low_end))/(high_end - low_end));

% Bin the data
output = discretize(scaled_raitings,[0 25 50 75 100]); % here are the bin limits

end