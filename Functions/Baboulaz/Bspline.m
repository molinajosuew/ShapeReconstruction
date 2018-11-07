function [b_n_m] = Bspline(n,m,PlotOrNot)

% Generate a discrete B-spline of order n and scale m
%
% Example: [b_n_m] = Bspline(3, 16,1)
%
% Parameters:
%   n               order of the B-spline
%   m               expansion factor (m): Resolution of the desired spline
%                   (# of taps in ]0 1]).
%   PlotOrNot       plot intermediate results or graphs, or not
%
% Output:
%   b_n_m           Discrete B-spline

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

% Generate the box spline at the desired scale
b_zero_m = [];
for i=1:m,
    b_zero_m = [b_zero_m 1];
    if n==0,
        b_n_m = b_zero_m;
        b_n_m = [b_n_m 0];
    end
end

if n>0,

% Select the transfer function of the desired B-spline
    [b_n_1 c_n_1] = Bspline_TransferFunction(n);

% First case: n odd ; m even
    if mod(m,2)==0,
        if mod(n,2)==1,
%             fprintf('\nn odd ; m even\n')
            temp = conv(b_zero_m,b_zero_m);
            for i=1:n-1,
                temp = conv(temp,b_zero_m);
            end

            b_n_m = conv(1/m^n*temp,b_n_1);

% Second case: n even ; m even
        else
%             fprintf('\nn even ; m even\n')
            temp = conv(b_zero_m,b_zero_m);
            for i=1:n-1,
                temp = conv(temp,b_zero_m);
            end

            b_n_m = conv(1/m^n*temp,c_n_1);
        end

        b_n_m = [0 b_n_m 0];
% Third case: n odd/even ; m odd
    else
%         fprintf('\nn odd/even ; m odd\n')
        temp = conv(b_zero_m,b_zero_m);
        for i=1:n-1,
            temp = conv(temp,b_zero_m);
        end

        b_n_m =  conv(1/m^n*temp,b_n_1);
        
        if mod(n,2)==0,
            b_n_m = [b_n_m 0];
        else
            b_n_m = [0 b_n_m 0];
        end
    end
end

if PlotOrNot
    support = linspace(-(n+1)/2,(n+1)/2,length(b_n_m));
    figure
%     stem_handles = stem(support,b_n_m,'LineWidth',1);
%     hold on
    plot_handles = plot(support,b_n_m,'LineWidth',2);
    hold off
    title(sprintf('B-Spline of order %g and scale %g',n,m),'Fontsize',16);
    set(gca,'Fontsize',16);
end  



