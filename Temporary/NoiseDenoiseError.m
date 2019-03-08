clear;
clc;

for j = 25 : - 5 : 25
    error_n = [];
    error_d = [];
    for i = 1 : 20
        load("./Images/Temporary/img" + i + ".mat");
        J = awgn(ConvolveAndDownsample(I, 2 * 4 - 1, 1, 1), j);
        K = abs(I - NonSepBernImg(NonSepAlg(J, 4), 4, size(I, 1), size(I, 2)));
        error_n = [error_n, sum(sum(K)) / numel(K)];
        K = abs(I - NonSepBernImg(NonSepAlg(wiener2(J), 4), 4, size(I, 1), size(I, 2)));
        error_d = [error_d, sum(sum(K)) / numel(K)];
    end
    save("./Statistics/error_" + j + ".mat", "error_n");
end