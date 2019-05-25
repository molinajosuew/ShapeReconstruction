clear;
clc;

load('/home/wjm/Documents/GitHub/ShapeReconstruction/Data/PosterImages/02/Matlab.mat');

figure;
imshow(I);
figure;
imshow(D);
figure;
imshow(R);

figure;
imagesc(I);
axis off;
figure;
imagesc(D);
axis off;
figure;
imagesc(R);
axis off;