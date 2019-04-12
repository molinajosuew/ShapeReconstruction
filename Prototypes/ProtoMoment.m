function M = ProtoMoment(K, i, j, J)
    if 0 <= i && 0 <= j
        KC = (size(K, 2) + 1) / 2;
        C1 = K(1 + j, KC - (size(J, 1) - 1) / 2 : KC + (size(J, 1) - 1) / 2);
        C2 = K(1 + i, KC - (size(J, 2) - 1) / 2 : KC + (size(J, 2) - 1) / 2);
        M = sum(sum(C1' * C2 .* J));
    else
        M = 0;
    end
end