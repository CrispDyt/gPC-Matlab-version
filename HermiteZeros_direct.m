function z = HermiteZeros(degree)
%
% HermiteZeros.m - Evaluates the zeros of Hermite polynomial
%                  with degree 'degree'.
%                  This routine does not employ the iterative solver. Instead,
%                  it takes advantage of the roots solve of matlab.
%                  The algorithm is much faster and stable for high order 
%                  Hermite polynomials.
%
% Syntax :    z = HermiteZeros( degree )
%
% Input  :    degree - degree of Hermite polynomial (number of zeros)
%
% Output :    z - zeros as array (degree x 1)
%
% By Dongbin Xiu   12/09/2002
%

coef = HermiteCoef( degree );

z = roots(fliplr(coef'));
