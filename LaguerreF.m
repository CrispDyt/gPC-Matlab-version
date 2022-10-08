function poly = LaguerreF(x, n, alpha)
%
% LaguerreF.m - Compute the Laguerre polynomial L^n(x;alpha)
%
% Syntax:   poly = LaguerreF(x, n, degree);
%
% Input :   x = x-coordinate in matrix form
%           n = order of polynomial
%           alpha > -1 (not mandatory) is the parameter.
%
% Output:   poly = polynomial values at location x stored in same form as x.
%
% by Dongbin Xiu   5/02/2003.
%

switch n
  case 0
    poly = ones(size(x));
  otherwise
    poly = polyval(LaguerreCoef(n,alpha), x);
end
