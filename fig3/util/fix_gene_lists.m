function plasmid_info = fix_gene_lists(plasmid_info,columns)

for i = 1:length(columns)
    temp = cellfun(@(x) strrep(x,']',''),plasmid_info{:,columns{i}},'UniformOutput',false);
    temp = cellfun(@(x) strrep(x,'[',''),temp,'UniformOutput',false);
    temp = cellfun(@(x) strrep(x,'''',''),temp,'UniformOutput',false);
    temp = cellfun(@(x) strrep(x,' ',''),temp,'UniformOutput',false);
    plasmid_info{:,columns{i}} = cellfun(@(x) strsplit(x,','),temp,'UniformOutput',false);
end

end

