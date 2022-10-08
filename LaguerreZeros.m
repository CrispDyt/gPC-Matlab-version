function z = LaguerreZeros(n, alpha)
%
% LaguerreZeros.m - Evaluates the zeros of Laguerre polynomial L^n(x;\alpha).
%                  This routine does not employ the iterative solver. Instead,
%                  it takes advantage of the roots solve of matlab.
%
% Syntax :    z = LaguerreZeros( n, alpha )
%
% Input  :    n - degree of Laguerre polynomial (number of zeros)
%             alpha > -1 is the parameter (this is not mandatory).
%
% Output :    z - zeros as array (degree x 1)
%
% By Dongbin Xiu   5/02/2003
%

z = roots(LaguerreCoef(n,alpha));
