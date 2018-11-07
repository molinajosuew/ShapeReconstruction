function [Dual,s_k] = DualSpline(n, m,MaxSupport,PlotOrNot)


% Calculate the discrete dual B-spline of order n, scale m and with maximal
% support MaxSupport.
%
% Example: [Dual,s_k] = DualSpline(3, 16,20,1)
%
% Parameters:
%   n               order of the dual B-spline
%   m               expansion factor (m): Resolution of the desired spline
%                   (# of taps in ]0 1]).
%   MaxSupport      support on which the infinite dual B-spline is computed
%                   (support centered in zero)
%   PlotOrNot       plot intermediate results or graphs, or not
%
% Output:
%   Dual            Discrete dual B-spline
%   s_k             impulse response of the filter for the Dual Spline

% --------------------------------------------------------------------- %
% Matlab code and data to reproduce results from paper                  %
% "Exact Feature Extraction using Finite Rate of Innovation Principles  %
% with an Application to Image Super-resolution"                		%
% by L. Baboulaz and P.L. Dragotti.                                     %
% Available at: http://www.commsp.ee.ic.ac.uk/~pld/                     %
%                                                                       %
% Author: Lo�c Baboulaz                                                 %
% Copyright (C) 2007 L. Baboulaz and P.L. Dragotti                      %
% Communications & Signal Processing Group,                             %
% Imperial College London, Exhibition Road, United Kingdom              %
%                                                                       %
% This program is free software; you can redistribute it and/or modify  %
% it under the terms of the GNU General Public License as published by  %
% the Free Software Foundation; either version 2 of the License, or (at %
% your option) any later version. This software is distributed in the   %
% hope that it will be useful, but without any warranty; without even   %
% the implied warranty of merchantability or fitness for a particular   %
% purpose.                                                              %
% See the GNU General Public License for more details                   %
% (enclosed in the file GPL).                                           %
%                                                                       %
% Latest modifications: 15 January 2008                                 %
% --------------------------------------------------------------------- %


if n==0,
    Dual = Bspline(n,m,PlotOrNot);
    support = -MaxSupport:1/m:MaxSupport;
    Dual = [zeros(1,(length(support)-length(Dual))/2) Dual zeros(1,(length(support)-length(Dual))/2)];
    s_k = [];
else
    % Impulse response support of the desired Dual Spline
    k = -MaxSupport:1:MaxSupport;
    support = k(1):1/m:k(end);

    % Generate the corresponding B-spline of order n
    [b_n_m] = Bspline(n,m,0);

    % Support of the B-spline
    BsplineSupport = n+1;

    % Find the transfer function coefficients of the filters for B-spline of
    % order 2n+1
    [b_n_1 c_n_1] = Bspline_TransferFunction(2*n+1);

    % Find the common denominator of the coefficients of the spline of order
    % 2n+1
    denom = 1/b_n_1(1);

    % Find the poles of the transfer function of the filter for Dual spline of
    % order n
    poles = roots(b_n_1*denom);
    poles = sort(poles,'descend');
    poles = poles(1:n);

    % Extend the support necessary to achieve the desired level of
    % precision
    Precision = 10^-16;
    MaxPole = max(poles);
    kprime = MaxSupport;
    while abs(MaxPole^kprime)>Precision,
        kprime = kprime+1;
    end
    kprime = -kprime:1:kprime;

    % Find the impulse response s_k of the filter for the Dual Spline
    weight = 1;
    ConvExp = [1];
    for i=1:n,
        weight = weight*poles(i)/(1-poles(i)^2);
        ConvExp = conv(ConvExp,poles(i).^abs(kprime));
    end

    % Recenter the support to the desired one
    IndexZero = (length(ConvExp)+1)/2;
    ConvExp = ConvExp(IndexZero-MaxSupport:IndexZero+MaxSupport);
    s_k = -denom*weight*ConvExp;
    if mod(n,2)==0, s_k = -s_k;end

    % Find Dual Spline from the linear combination of B-splines
    Dual = zeros(1,m*(length(k)-1 ) +1);
    Initial_Bspline = [b_n_m zeros(1,m*(length(k)-1 -BsplineSupport) )];

    % BsplineSupport
    for i=1:length(k)-BsplineSupport,
        Dual = Dual + s_k(round(i+BsplineSupport/2))*circshift(Initial_Bspline, [0 (i-1)*m]);
    end
    if mod(n,2)==0, Dual = [0 circshift(Dual(2:end-1),[0 round(0.5*m)]) 0];end
end


if PlotOrNot
    figure;
%     stem_handles = stem(support,Dual,'LineWidth',1);
%     hold on
    plot_handles = plot(support,Dual,'LineWidth',2);
    hold off
    title(sprintf('Dual spline of order %g and scale %g',n,m),'Fontsize',16);
    set(gca,'Fontsize',16);
end

