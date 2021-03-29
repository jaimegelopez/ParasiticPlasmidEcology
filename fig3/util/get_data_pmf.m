function [pmf,xvec,indvec] = get_data_pmf(obs_vec,count_cutoff)

%This function gets an empirical pmf from a dataset

[counts] = count_elements(obs_vec,max(obs_vec));

total_count = sum(counts);

if count_cutoff > 0
    below_cutoff = find(counts < count_cutoff);
    below_cutoff = below_cutoff(1) -1;
    counts = counts(1:below_cutoff);
end

pmf = counts/total_count;

xvec = 0:(length(pmf)-1);
indvec = 1:length(pmf);

end

