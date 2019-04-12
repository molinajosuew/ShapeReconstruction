function M = ProtoMoment(C_x, C_y, i, j, J)
    if 0 <= i && 0 <= j
        M = (i + r) * sum(sum(C(i + r + 0, C_y)' * C(j + s + 1, C_x) .* J2));
        KC = (size(C, 2) + 1) / 2;
        C_y = C(1 + j, KC - (size(J, 1) - 1) / 2 : KC + (size(J, 1) - 1) / 2);
        C_x = C(1 + i, KC - (size(J, 2) - 1) / 2 : KC + (size(J, 2) - 1) / 2);
        M = sum(sum(C_y' * C_x .* J));
    else
        M = 0;
    end
end