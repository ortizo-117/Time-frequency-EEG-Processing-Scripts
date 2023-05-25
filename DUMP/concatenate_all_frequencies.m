function output = concatenate_all_frequencies(gamma_in,lf_in)

%% for testing
my_names = fieldnames(gamma_in);

for i = 1:length(my_names)-2
    my_gamma = eval(['gamma_in.' char(my_names(i)) '(2:end,:);']);
    my_lf = eval(['lf_in.' char(my_names(i)) ';']);
    data_in = [my_lf;my_gamma];
    
    eval(['output.' char(my_names(i)) ' = data_in;']);
end

output.x = gamma_in.x;
output.y = [lf_in.y(1:end-1) gamma_in.y];



end
