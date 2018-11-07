function [SF_coef] = StrangFix_coefficients(SplineOrder,Scale,Support,PlotOrNot)

% Calculate the coefficients necessary for polynomial reproduction with
% B-Spline. These coefficients are referred as the Strang-Fix coefficients
%
% Example: [SF_coef] = StrangFix_coefficients(3,8,[-20 20],512,1)
%
% Parameters:
%   SplineOrder:    order of the B-Spline used for polynomial reproduction
%   Scale:          Scale of the discrete B-Spline (#pixel per unit)
%   Support:        1x2 matrix of the support considered for polynomial
%                   reproduction
%   PlotOrNot       plot intermediate results or graphs
%
% Output:
%   SF_coef         the Strang-Fix coefficients for polynomial reproduction

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


% Support of Dual spline is 10 times the support for reproduction
DualSupport = round(10*Support(2));
Dual = DualSpline(SplineOrder,Scale,DualSupport,PlotOrNot);


% Corresponding time frame
Time = Support(1)-DualSupport:1/Scale:Support(2)+DualSupport+(Scale-1)/Scale;


% Compute the Strang-Fix coefficients
for Degree=0:SplineOrder,
    % Polynomial to approximate
    Polynomial = (Time).^Degree;
    % Compute the inner product between the signal and the Dual spline
    Coef=[];
    for k=0:Support(2)-Support(1),
        Coef(k+1,1) = Polynomial(1 ,k*Scale+1:k*Scale+size(Dual,2))...
        *Dual'*(1/Scale);
    end
    SF_coef(Degree+1,:) = Coef';
    
    if PlotOrNot,
        % Reconstruct the polynomial from the coefficients and the spline only.
        TimeRepro = Support(1)-(SplineOrder+1)/2 :1/Scale :Support(2)+(SplineOrder+1)/2; %-1/Scale
        PolyRepro = zeros(size(TimeRepro,2),1);
        Spline = Bspline(SplineOrder,Scale,0);

        % plot the scaled shifted versions of the B-Spline
        figure;
        for k=0:Support(2)-Support(1),
            plot(TimeRepro(1,k*Scale+1:k*Scale + size(Spline,2))',Coef(k+1,1)*Spline', 'b','LineWidth',2);
            hold on;
            PolyRepro(k*Scale+1:k*Scale + size(Spline,2),1) = ...
                PolyRepro(k*Scale+1:k*Scale + size(Spline,2),1) + ...
                Coef(k+1,1)*Spline';
        end
        % Plot the reproduction and the original polynomial
        plot(TimeRepro,TimeRepro.^Degree,':k','LineWidth',2)
        hold on
        plot(TimeRepro,PolyRepro,'r','LineWidth',2);
        title(sprintf('Reproduction of t^{%g} with shifted B-splines of order %g',Degree,SplineOrder),'Fontsize',16);
        set(gca,'Fontsize',16);
        axis('tight');
    end
end

