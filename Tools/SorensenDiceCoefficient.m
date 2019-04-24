function DSC = SorensenDiceCoefficient(I, R)
    TP = sum(sum((I == R) .* I));
    FP = sum(sum((I ~= R) .* R));
    FN = sum(sum((I ~= R) .* I));
    
    DSC = 2 * TP / (2 * TP + FP + FN);
end