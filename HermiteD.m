function polyd = HermiteD(x, degree)
%
% HermiteD.m - Compute the derivative of Hermite polynomial H^n(x)
%
% Syntax:   polyd = HermiteD(x, degree);
%
% Input :   x = x-coordinate in matrix form
%           degree = order of polynomial
%
% Output:   polyd = derivative values at location x stored in same form as x.
%
% by Dongbin Xiu   5/03/2002
%
switch degree
  case 0
    polyd = zeros(size(x));
  case 1
    polyd = ones(size(x));
  otherwise
    polyd = degree*HermiteF(x, degree-1);
end
