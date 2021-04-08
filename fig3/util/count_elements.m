function [counts,elements] = count_elements(vec,max_plasmids)

%Counts occurences of different elements within a sample

elements = 0:max_plasmids;
counts = zeros(size(elements));

for i = 1:length(elements)
    counts(i) = sum(vec == elements(i));
end


end

