% Input
% x : convolution function : 1 by n matrix with n odd (TEMPORARY ASSUMPTION)
% y : convolution function : 1 by m matrix
% s : step size : positive integer

% Output
% z : convolution : 1 by k matrix

function z = ProtoToeplitzConvolution1(x, y, s, f)
    c = (size(x, 2) + 1) / 2;
    if f
        if floor(c) ~= c
            disp("convolution function has no center");
            return;
        end
    else
        while c + s <= size(x, 2)
            c = c + s;
        end
    end
    i = 0;
    r = [zeros(1, min(- min(c - i * s, 1) + 1, size(y, 2))), x(max(c - i * s, 1) : c - i * s + min(size(x, 2) - (c - i * s), size(y, 2) - 1)), zeros(1, size(y, 2) - min(size(x, 2) - (c - i * s), size(y, 2) - 1) - 1)];
    while nnz(r)
        T(i + 1, :) = r;
        i = i + 1;
        r = [zeros(1, min(- min(c - i * s, 1) + 1, size(y, 2))), x(max(c - i * s, 1) : c - i * s + min(size(x, 2) - (c - i * s), size(y, 2) - 1)), zeros(1, size(y, 2) - min(size(x, 2) - (c - i * s), size(y, 2) - 1) - 1)];
    end
    z = transpose(T * transpose(y));
    T
    transpose(y)
end