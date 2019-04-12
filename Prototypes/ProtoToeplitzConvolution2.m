% Input
% x : convolution function : 1 by n matrix
% y : convolution function : 1 by m matrix
% s : step size : positive integer
% is_valid : true for a valid convolution : boolean

% Output
% z : convolution result : 1 by k matrix

function z = ProtoToeplitzConvolution2(x, y, step, type)
    t = [x, zeros(1, step * size(y, 2) - step)];
    T = zeros(size(t, 2), size(y, 2));
    for i = 1 : size(y, 2)
        T(:, i) = circshift(t, step * (i - 1));
    end
    if nargin == 4 && type ~= "normal"
        if type == "valid"
            T = T(find(T(:, 1) == 0, 1) - 1 : find(T(:, end) ~= 0, 1), :);
        elseif type == "centered"
            c = (size(x, 2) + 1) / 2;
            if floor(c) ~= c
                disp("function with no center given");
                return;
            end
            T = T(find(T(:, 1) == x(c), 1) : find(T(:, end) == x(c), 1), :);
        else
            disp("invalid type");
            return;
        end
    end
    z = transpose(T * transpose(y));
end