clear;
clc;

n = 4;

ra = - 2;
rb = 6;
ca = - 2;
cb = 6;

rQ = 1;
cQ = 1;

rN = rQ * (rb - ra) + 1;
cN = cQ * (cb - ca) + 1;

I = RandomStablyBounded4thDegreeBivariatePolynomial(ra, rb, ca, cb, rN, cN, 250000);
figure;
imshow(I);

[C, D] = ProtoPowerReconstruction(I, n, ra, rb, ca, cb);

R = RealizePolynomial(C, rb, ra, ca, cb, rN, cN);
figure;
imshow(D);
figure;
imshow(R);

% I was testing small images to see if the C's align.