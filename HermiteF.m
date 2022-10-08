function poly = HermiteF(x, degree)
%
% HermiteF.m - Compute the Hermite polynomial H^n(x)
%
% Syntax:   poly = HermiteF(x, degree);
%
% Input :   x = x-coordinate in matrix form
%           degree = order of polynomial
%
% Output:   poly = polynomial values at location x stored in same form as x.
%
% by Dongbin Xiu   5/03/2002
%

switch degree
  case 0
    poly = ones(size(x));
  case 1
    poly = x;
  otherwise
    poly = x.*HermiteF(x, degree-1) - (degree-1)*HermiteF(x, degree-2);
end
