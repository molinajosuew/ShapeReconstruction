% I - Image
% n - B-Spline Degree
% m - B-Spline Dilation Factor
% t - Downsampling Factor

function ConvolveAndDownsampleOut = ConvolveAndDownsample(I, n, m, t)
    helper = conv2(wjmBSpline(n, m), wjmBSpline(n, m), I, 'same');

%     % This is temporarily done due to not being able to down-sample in the convolution above.
%     ConvolveAndDownsampleOut = imresize(helper, 1 / t);
    
    ConvolveAndDownsampleOut = helper(1 : t : end, 1 : t : end);

%     ConvolveAndDownsampleOut = zeros(size(I, 1) / t, size(I, 2) / t);
%     for i = 1 : size(I, 1) / t
%         for j = 1 : size(I, 2) / t
%             ConvolveAndDownsampleOut(i, j) = helper(i * t, j * t);
%         end
%     end
end