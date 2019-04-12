% I - Image
% n - B-Spline Degree
% m - B-Spline Dilation Factor
% t - Downsampling Factor

function ConvolveAndDownsampleOut = ConvolveAndDownsample(I, n, m, t)
    helper = conv2(wjmBSpline(n, m), wjmBSpline(n, m), I, 'same');
    ConvolveAndDownsampleOut = helper(1 : t : end, 1 : t : end);
end