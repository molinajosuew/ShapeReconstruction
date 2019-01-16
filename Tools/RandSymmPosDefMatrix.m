function M = RandSymmPosDefMatrix(n)
    M = (2 * rand(n) - 1) * 100;
    for i = 1 : n
        M(i, i) = abs(M(i, i));
    end
    M = M * M';
end