function plasmid_info = assign_taxonomy(plasmid_info)

%Assigns taxonomy to an organism

organism = plasmid_info.organism;

splits = cellfun(@(x) strsplit(x,' '), organism,'UniformOutput',false);

genus = cellfun(@(x) lower(x{1}), splits,'UniformOutput',false);
species = cellfun(@(x) lower(join(x(1:2),' ')), splits,'UniformOutput',false);
species = cellfun(@(x) x{1}, species,'UniformOutput',false);

plasmid_info.genus = genus;
plasmid_info.species = species;

end

