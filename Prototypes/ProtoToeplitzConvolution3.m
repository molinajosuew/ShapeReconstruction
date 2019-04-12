function z = ProtoToeplitzConvolution3(x, y, s)
    T(:, 1) = [x, zeros(1, size(y, 2) * s - s)]';
    for i = 2 : size(y, 2)
        T(:, i) = circshift(T(:, i - 1), s);
    end
    z = transpose(T * transpose(y));
end