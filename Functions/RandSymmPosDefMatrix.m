function M = RandSymmPosDefMatrix(n)
    M = 2 * rand(n) - 1;
    for i = 1 : n
        M(i, i) = abs(M(i, i));
    end
    M = M * M';
end