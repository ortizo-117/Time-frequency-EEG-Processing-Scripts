function output = concatenate_all_results(gamma_results,lf_results)

%% for testing
my_names = fieldnames(gamma_results);

for i = 1:length(my_names)
    my_gamma = eval(['gamma_results.' char(my_names(i)) ';']);
    my_lf = eval(['lf_results.' char(my_names(i)) ';']);
    data_in = [my_lf;my_gamma];
    
    eval(['output.' char(my_names(i)) ' = data_in;']);
end

end


