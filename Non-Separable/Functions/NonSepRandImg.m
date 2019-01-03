function [NonSepRandImgOut, Coeffs] = NonSepRandImg(n, a, b)
NonSepRandImgOut = zeros(a, b);
while NonSepRandImgCheck(NonSepRandImgOut)
    NonSepRandImgOut = zeros(a, b);
    Coeffs = [];
    for i = 0 : n
        for j = 0 : n - i
            coeff = (2 * rand() - 1) * 100;
            Coeffs = union(Coeffs, coeff);
            NonSepRandImgOut = NonSepRandImgOut + coeff * NonSepBern(1, a, b, n, i, j);
        end
    end
    NonSepRandImgOut = im2double(NonSepRandImgOut <= 0);
end
end