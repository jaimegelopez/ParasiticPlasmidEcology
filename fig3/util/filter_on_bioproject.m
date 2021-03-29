function filtered_origin_table = filter_on_bioproject(origin_table)

%This function filters genomes based on the description of their bioproject

bioproject_data = readtable('bioproject_list_fixed.csv','ReadRowNames',true);

or_bioprojects = origin_table.bioproject;
valid_bioprojects = intersect(or_bioprojects,bioproject_data.Properties.RowNames);

absent_ids = setdiff(unique(or_bioprojects),valid_bioprojects);
absent_ind = false(size(origin_table,1),1);
for i = 1:length(absent_ids)
    absent_ind = absent_ind | strcmp(origin_table.bioproject,absent_ids{1});
end

sub_bioproj_info = bioproject_data(valid_bioprojects,:);

sus_text = {'engineer','clone','cloni','artificial','laboratory','recomb',...
    'industr','chimera','Evolution of bacterial fitness with an expanded genetic code'};

sus_bioproj = sub_bioproj_info(contains(sub_bioproj_info.ProjectTitle,sus_text,'IgnoreCase',true),:);

sus_index = absent_ind;
for i = 1:size(sus_bioproj,1)
    sus_id = sus_bioproj.Properties.RowNames{i};
    
    sus_index = sus_index | strcmp(origin_table.bioproject,sus_id);

end
    
filtered_origin_table = origin_table(~sus_index,:);

filtered_origin_table.title = bioproject_data{filtered_origin_table.bioproject,'ProjectTitle'};

end

