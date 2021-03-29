function reordered_combos =  reorder_combos(combos)

reordered_combos = zeros(size(combos));

for i = 1:size(combos,1)
    combo = combos(i,:);
    reordered_combos(i,:) = [min(combo),max(combo)];
end

end

