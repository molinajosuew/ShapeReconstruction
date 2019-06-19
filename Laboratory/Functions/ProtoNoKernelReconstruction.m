function C = ProtoNoKernelReconstruction(I, n)
    M = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
    row = 1;
    for r = 0 : ceil(n / 2)
        for s = 0 : ceil(n / 2)
            col = 1;
            for k = 0 : n
                for l = 0 : n - k
                    M(row + 0, col) = (k + r) * ProtoNoKernelMoment(I, k + r - 1, l + s);
                    M(row + 1, col) = (l + s) * ProtoNoKernelMoment(I, k + r, l + s - 1);
                    col = col + 1;
                end
            end
            row = row + 2;
        end
    end
    C = quadprog(M' * M, [], [], [], eye(1, size(M, 2)), 1, [], [], [], optimset('display', 'off'))';
end