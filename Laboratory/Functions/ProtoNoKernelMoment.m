function M = ProtoNoKernelMoment(I, i, j)
    if i >= 0 && j >= 0
        [X, Y] = meshgrid(0 : size(I, 1) - 1, 0 : size(I, 2) - 1);
        M = sum(sum(X .^ i .* Y .^ j .* I));
    else
        M = 0;
    end
end