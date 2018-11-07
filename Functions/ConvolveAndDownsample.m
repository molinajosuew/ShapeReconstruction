% I - Image
% n - B-Spline Degree
% m - B-Spline Dilation Factor
% t - Downsampling Factor

function ConvolveAndDownsampleOut = ConvolveAndDownsample(I, n, m, t)
% ConvolutionAndDownsamplingOut = imresize(conv2(BSpline(n, m), BSpline(n, m), I, 'same'), 1 / t, 'box');
helper = conv2(BSpline(n, m), BSpline(n, m), I, 'same');
ConvolveAndDownsampleOut = zeros(size(I, 1) / t, size(I, 2) / t);
for i = 1 : size(I, 1) / t
    for j = 1 : size(I, 2) / t
        ConvolveAndDownsampleOut(i, j) = helper(i * t, j * t);
    end
end
end