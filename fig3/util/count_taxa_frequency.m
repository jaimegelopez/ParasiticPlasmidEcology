function [taxa,counts] = count_taxa_frequency(tbl,level)

%Counts frequency of the different taxa

taxa = unique(tbl{:,level});
counts = zeros(size(taxa));
for i = 1:length(taxa)
    counts(i) = size(tbl(strcmp(tbl{:,level},taxa{i}),:),1);
end

[counts,idx] = sort(counts,'descend');
taxa = taxa(idx);

end