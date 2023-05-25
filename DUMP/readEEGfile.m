% Author: Oscar Ortiz
% Date: 2021/05/21
% Purpose:Read .csv containing EEGdata, and extract data and channel names
% Inputs: filename - Name of the .csv file corresponding to the EEG data
% Outputs: data_out = matrix containing EEG time-series. Each column
%   contains one channel
%          chan_names = Cell vector containing the name of each channel.


function [data_out, chan_names] = readEEGfile(filename)


T = readtable(filename);

data_out = T{:,:};

chan_names = T.Properties.VariableNames;


end
