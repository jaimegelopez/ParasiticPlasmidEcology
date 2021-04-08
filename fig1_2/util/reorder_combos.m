function reordered_combos =  reorder_combos(combos)

%This script reorders combinations to always be in the same ascending
%order

reordered_combos = zeros(size(combos));

for i = 1:size(combos,1)
    combo = combos(i,:);
    reordered_combos(i,:) = [min(combo),max(combo)];
end

end

