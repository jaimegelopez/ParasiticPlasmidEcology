function origin_table = get_origin_table(plasmid_info)

%This function takes a list of elements and aggregates it by origin

%Construct the table
origin_list = unique(plasmid_info.Origin);
column_list = {'n_plasmids','n_chromosomes','n_other','cas_present'};
temp = nan(length(origin_list),length(column_list));
origin_table = array2table(temp,'RowNames',origin_list,...
    'VariableNames',column_list);
origin_table.organism = cell(size(origin_table.n_plasmids));

%Loop through origins
for i = 1:length(origin_list)
    
    %Select origin i and get its elements
    origin = origin_table.Properties.RowNames{i}; 
    subtable = plasmid_info(strcmp(plasmid_info.Origin,origin),:);
    
    %Count the number of plasmids, chromosomes, and other element classes
    origin_table.n_plasmids(i) = ...
        sum(contains(subtable.Assigned_Molecule_Location_Type,'Plasmid'));
    origin_table.n_chromosomes(i) = ...
        sum(contains(subtable.Assigned_Molecule_Location_Type,'Chromosome'));
    origin_table.n_other(i) = size(subtable,1) - ...
        origin_table{origin,'n_chromosomes'} - origin_table{origin,'n_plasmids'};
    
    origin_table.organism{i} = subtable.Organism{1};
    origin_table.bioproject{i} = subtable.Bioproject{1};

    
    %Count cas gene presence
    chrom_subtable = ...
        subtable(strcmp(subtable.Assigned_Molecule_Location_Type,'Chromosome'),:);
    origin_table.cas_present(i) = ...
        sum(~cellfun(@(x) isempty(x{1}),subtable.cas))>0;
    origin_table.cas_count(i) = ...
        sum(cellfun(@(x) sum(~strcmp(x,'')),subtable.cas));
    origin_table.cas_present_chrom(i) = ...
        sum(~cellfun(@(x) isempty(x{1}),chrom_subtable.cas))>0;
end

