function M = ProtoNonSeparableBernsteinMoment(n, i, j, J, C, C_x, C_y)
    if n >= 0 && i >= 0 && j >= 0
        M = sum(sum(reshape(C(1 + n, 1 + i, 1 + j, C_y, C_x), size(C_y, 2), size(C_x, 2)) .* J));
    else
        M = 0;
    end
end