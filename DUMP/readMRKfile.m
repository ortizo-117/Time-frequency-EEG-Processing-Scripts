% Author: Oscar Ortiz
% Date: 2021/05/21
% Purpose:Read MRK containing the raitings, and extracts the recorded
%   raitings and exports ne number of raitings recorded as well as a
%   numberical vector containing the recorded ratings 
% Inputs: filename - Name of the MRK file corresponding to the experimental
%   condition.
% Outputs: extracted_raitings = vector of extracted raitings.
%          n_raitings = number of raitings.


function [extracted_raitings, n_raitings] = readMRKfile(filename)

text = fileread(filename);

expr = '[^\n]*Comment[^\n]*';

matches = regexp(text,expr,'match');

extracted_raitings = zeros(1,length(matches));


for i = 1:length(matches)
    cur_stim = char((matches(i)));
    temp = split(cur_stim,"t");
    temp2 = temp(2);
    temp3 = char(temp2);
    temp4 = regexprep(temp3, '\s+', '');
    cur_rat = temp4(1);
    extracted_raitings(i) = str2double(cur_rat);
end

n_raitings = length(extracted_raitings);

end
