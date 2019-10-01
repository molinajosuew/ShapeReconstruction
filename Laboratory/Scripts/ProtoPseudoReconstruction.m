clear;
clc;

I = GetImageOfPower(GetRandomPower(4, - 1, 1, - 1, 1, false), - 1, 1, - 1, 1, 2 ^ 9 + 1, 2 ^ 9 + 1);

P = [6435/8388608, -7293/8388608, -8415/1048576, 9945/1048576, 85085/2097152, -109395/2097152, -153153/1048576, 255255/1048576, 3828825/4194304, 3828825/4194304, 255255/1048576, -153153/1048576, -109395/2097152, 85085/2097152, 9945/1048576, -8415/1048576, -7293/8388608, 6435/8388608];
n = 4;
B = ScalingIntegers(P);
D = im2double(conv2(B, B, I));

for k = 0 : n + ceil(n / 2)
    K(1 + k, :) = PseudosplinesPolCoeffs(n, k, P, - 1, size(I, 1)); % not optimized
end

M = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
row = 1;

for r = 0 : ceil(n / 2)
    for s = 0 : ceil(n / 2)
        col = 1;
        
        for i = 0 : n
            for j = 0 : n - i
                if 0 <= i + r - 1
                    M(row + 0, col) = (i + r) * sum(sum(K(j + s + 1, 1 : size(D, 1))' * K(i + r + 0, 1 : size(D, 1)) .* D));
                else
                    M(row + 0, col) = 0;
                end
                
                if 0 <= j + s - 1
                    M(row + 1, col) = (j + s) * sum(sum(K(j + s + 0, 1 : size(D, 1))' * K(i + r + 1, 1 : size(D, 1)) .* D));
                else
                    M(row + 1, col) = 0;
                end
                
                col = col + 1;
            end
        end
        
        row = row + 2;
    end
end

R = GetImageOfPower(quadprog(M' * M, [], [], [], eye(1, size(M, 2)), 1, [], [], [], optimset('display', 'off'))', 0, size(I, 1), 0, size(I, 2), 2 ^ 9 + 1, 2 ^ 9 + 1);
imshow(I);
figure;
imshow(abs(I - R));