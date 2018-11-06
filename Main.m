J = GetSeparableBernsteinImage(Algorithm(ConvolveAndDownsample(I, 8, 1, 1), 4), 4, 4, size(I, 1), size(I, 2));
figure;
imshow(I);
figure;
imshow(J);
figure;
imshow(abs(I - J));