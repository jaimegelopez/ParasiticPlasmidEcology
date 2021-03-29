function [N,K,types,type_combo_index] = generate_plasmid_matrix(n_plasmid)

%This function generates the matrices for the multi-plasmid model

%Generate all possible cell types
types = dec2bin(2^n_plasmid-1:-1:0)-'0'; 

%Generate all non-redundant conjugation combos
type_combo_index = nchoosek(1:size(types,1),2);

%Get the number of types and combos
num_types = size(types,1);
num_combos = size(type_combo_index,1);

%Get each cell type involved in these combos
type_combos = cell(num_combos,2);
for i = 1:num_combos
    type_combos{i,1} = types(type_combo_index(i,1),:);
    type_combos{i,2} = types(type_combo_index(i,2),:);
end

manhat1 = @(v1,v2) sum(abs(v1 - v2)) == 1; 
more_plasmids = @(v1,v2) sum(v1) > sum(v2); 


%Initiate the conjugation matrix
N = zeros(num_types,num_combos);

%Loop through and define conjugation matrix
for i = 1:num_types
    %For each type, loop through combos and determine how each combo
    %affects that type (does it lose, gain, etc.)
    type = types(i,:);
    for j = 1:num_combos
        combo1 = type_combos{j,1};
        combo2 = type_combos{j,2};
        
        %Check if the type is involved in this conjugation 
        involved = isequal(type,combo1) | isequal(type,combo2);
        
        %If it is involved, figure out how many ways it can lose. One way
        %to lose for each plasmid it doesn't have that the other plasmid
        %does
        lost = 0;
        if involved && isequal(type,combo1)
            lost = sum((~type) & combo2);
        elseif involved && isequal(type,combo2)
            lost = sum((~type) & combo1);
        end
            
        %Check how many of the combo members are 1 away from the type and
        %also have fewer plasmids.
        one_less = [manhat1(type,combo1) & more_plasmids(type,combo1),...
            manhat1(type,combo2) & more_plasmids(type,combo2)];
        
        %Check, for the combos members are one away and have fewer plasmids,
        %if the other part of the combo has the missing plasmid
        produced = 0;
        for k = 1:length(one_less)
            if one_less(k) == 1
                
                %Get the sequences of the acceptor and the donor
                acceptor = type_combos{j,k};
                acceptor_index = zeros(size(one_less));
                acceptor_index(k) = 1;
                donor = type_combos{j,~acceptor_index};
                
                %Find the acceptor's missing plasmids
                mismatch = find((~acceptor) & type);
                
                %Increment produced if the donor has the missing plasmid
                produced = produced + donor(mismatch);
            end
        end
    
        N(i,j) = produced - lost;
        
    end
end

%Initialize loss matrix
K = zeros(num_types,num_types);

%Loop through and define loss matrix
for i = 1:num_types
    type = types(i,:);
    for j = 1:num_types
        other_type = types(j,:);
        if isequal(type,other_type)
            K(i,j) = -sum(type);
        elseif manhat1(type,other_type) && more_plasmids(other_type,type)
            K(i,j) = 1;
        end
    end
end


end

