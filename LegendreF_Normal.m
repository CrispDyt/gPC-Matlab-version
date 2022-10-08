function poly = LegendreF_Normal(x, n)
%
% LegendreF_Normal.m - Compute the normalized Legendre polynomial;
%
% Syntax:   poly = LegendreF_Normal(x, n);
%
% Input :   x = x-coordinate in matrix form
%           n = order of polynomial
%
% Output:   poly = polynomial values at location x stored in same form as x.
%
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% Algorithm: P_n(x) = sqrt(2n+1) L_n(x), where L_n(x) is the classical
%            Legendre polynomial.
%
% By Dongbin Xiu   06/13/2005
%

  poly = sqrt(2*n+1)*JacobiF(x,n,0,0);
