clear;
clc;

for i = 70 : - 5 : 30
    j = randi(20);
    load("/home/wjmolina/Documents/git-repos/ShapeReconstruction/Images/Temporary/img" + j + ".mat");
    J = awgn(ConvolveAndDownsample(I, 2 * 4 - 1, 1, 1), i);
    K = NonSepBernImg(NonSepAlg(J, 4), 4, size(J, 1), size(J, 2));
    L = NonSepBernImg(NonSepAlg(wiener2(J), 4), 4, size(J, 1), size(J, 2));
    save("/home/wjmolina/Documents/git-repos/ShapeReconstruction/Statistics/rec_" + j + "_" + i + ".mat", "I", "J", "K", "L");
end