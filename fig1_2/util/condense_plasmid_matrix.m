function [cN,cK,cond_combos,num_variants] = condense_plasmid_matrix(N,K,types,type_combo_index,n_plasmid)

%This function condenses the full multi-plasmid matrix into a form
%where all plasmids are identical

%Get number of variants for each type of plasmid containing cell, ranging
%from zero to n_plasmid plasmids
num_variants = arrayfun(@(x) nchoosek(n_plasmid,x),1:n_plasmid)';
num_variants = [1; num_variants];

%Get all combination of plasmid numbers
cond_combos =  nchoosek(0:n_plasmid,2);
cond_combos = [cond_combos; (1:n_plasmid)',(1:n_plasmid)'];
num_combos = size(cond_combos,1);

%Generate matrices
cN = zeros(n_plasmid + 1,num_combos);
cK = zeros(n_plasmid + 1, n_plasmid + 1);

%Get full combinations in terms of plasmid number
ctypes = sum(types,2);
ccombos = [ctypes(type_combo_index(:,1)),ctypes(type_combo_index(:,2))]; 

%Reorder all combos to be [lower,higher]
cond_combos = reorder_combos(cond_combos);
ccombos = reorder_combos(ccombos);

%Loop through all entries in the original conjugation matrix and use it to
%construct entries in the condensed matrix
tracker = 0;
for i = 1:length(ctypes)
    %Get index of condensed population
    pop_index = ctypes(i) + 1;
    for j = 1:size(ccombos,1)
        %Get index of condensed combo
        combo_index = find(cond_combos(:,1) == ccombos(j,1) & cond_combos(:,2) == ccombos(j,2));
        tracker = tracker + isempty(combo_index);
        cN(pop_index,combo_index) = cN(pop_index,combo_index) + N(i,j);    
    end
end

%Loops through original loss matrix and construct entries in condensed
%matrix
for i = 1:length(ctypes)
    %Get indexed of condensed population
    pop_index1 = ctypes(i) + 1;
    for j = 1:length(ctypes)
        %Get index of condensed combo
        pop_index2 = ctypes(j) + 1;
        cK(pop_index1,pop_index2) = cK(pop_index1,pop_index2) + K(i,j);    
    end
end

%Note: there are some combinatorial multipliers that are included
%downstream, this matrix works with the populations divided by their
%respective n choose k.


end

