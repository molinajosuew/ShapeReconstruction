function [b_n c_n] = Bspline_TransferFunction(n)

% Calculate the transfer function for generating discrete B-spline
%
% Example: [b_n c_n DesiredIndex_b_n DesiredIndex_c_n] = Bspline_TransferFunction(3)
%
% Parameters:
%   n               order of the B-spline
%
% Output:
%   b_n             Discrete B-spline of order n
%   c_n             Discrete B-spline of order n shifted by -1/2

% --------------------------------------------------------------------- %
% Matlab code and data to reproduce results from paper                  %
% "Feature Extraction for Image Super-resolution using Finite Rate of   %
% Innovation" by L. Baboulaz and P.L. Dragotti.                         %
% Available at: http://www.commsp.ee.ic.ac.uk/~pld/                     %
%                                                                       %
% Author: Loïc Baboulaz                                                 %
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


% B-spline centered in 0 of order n=1
b_1 = [1];
index_b_1 = [0];
% B-spline centered in -1/2 of order n=1
c_1 = [1/2 1/2];
index_c_1 = [-1 0];

% Indices k of the transfer functions
DesiredIndex_b_n = -fix(n/2):1:fix(n/2);
DesiredIndex_c_n = -fix((n-1)/2)-1:1:fix((n-1)/2);


% Number of zeros to add on each side of b_1
NbZeros = (length(DesiredIndex_b_n) - 1)/2;
b_1 = [zeros(1,NbZeros) b_1 zeros(1,NbZeros)];

% Number of zeros to add on each side of c_1
NbZeros = (length(DesiredIndex_c_n) - 2)/2;
c_1 = [zeros(1,NbZeros) c_1 zeros(1,NbZeros)];

if n==1,
    b_n = b_1;
    c_n = c_1;
else
    b_prev = b_1;
    c_prev = c_1;

    % Iterative approach to find the transfer functions
    for i=2:n,

        % B-spline centered in 0 of order n=1
        for k=1:length(DesiredIndex_b_n),
            if DesiredIndex_b_n(k)< DesiredIndex_c_n(1)||DesiredIndex_b_n(k)> DesiredIndex_c_n(end),
                cprev_k = 0;
            else
                ind = find(DesiredIndex_b_n(k)== DesiredIndex_c_n(:));
                cprev_k = c_prev(ind);
            end
            if DesiredIndex_b_n(k)-1< DesiredIndex_c_n(1)||DesiredIndex_b_n(k)-1> DesiredIndex_c_n(end),
                cprev_kMinusOne = 0;
            else
                ind = find(DesiredIndex_b_n(k)-1== DesiredIndex_c_n(:));
                cprev_kMinusOne = c_prev(ind);
            end

            b_i(k) = ((DesiredIndex_b_n(k) + (i+1)/2)*cprev_k + ((i+1)/2 - DesiredIndex_b_n(k))*cprev_kMinusOne) / i;
        end

        % B-spline centered in -1/2 of order n=1
        for k=1:length(DesiredIndex_c_n),
            if DesiredIndex_c_n(k)< DesiredIndex_b_n(1)||DesiredIndex_c_n(k)> DesiredIndex_b_n(end),
                bprev_k = 0;
            else
                ind = find(DesiredIndex_c_n(k)== DesiredIndex_b_n(:));
                bprev_k = b_prev(ind);
            end
            if DesiredIndex_c_n(k)+1< DesiredIndex_b_n(1)||DesiredIndex_c_n(k)+1> DesiredIndex_b_n(end),
                bprev_kPlusOne = 0;
            else
                ind = find(DesiredIndex_c_n(k)+1== DesiredIndex_b_n(:));
                bprev_kPlusOne = b_prev(ind);
            end

            c_i(k) = ((DesiredIndex_c_n(k) + (i+2)/2)*bprev_kPlusOne + (i/2 - DesiredIndex_c_n(k))*bprev_k) / i;
        end
        b_prev = b_i;
        c_prev = c_i;
    end

    b_n = b_i;
    c_n = c_i;
end