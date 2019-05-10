function NonSepRandImgCheckOut = NonSepRandImgCheck(I)
NonSepRandImgCheckOut = false;
if length(unique(union(union(union(I(1, :), I(end, :)), I(:, 1)), I(:, end)))) ~= 1 || length(unique(I)) == 1 || length(unique(triu(rot90(I)))) ~= 1
    NonSepRandImgCheckOut = true;
end
end