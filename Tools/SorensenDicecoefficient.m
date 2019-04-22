function DSC = SorensenDicecoefficient(X, Y)
    DSC = 2 * sum(sum(X == Y)) / (numel(X) + numel(Y));
end