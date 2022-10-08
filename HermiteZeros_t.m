function z = HermiteZeros_t(degree, t)
%
% HermiteZeros_t.m - Evaluates the zeros of Hermite polynomial h_n(x;t)
%                  with degree n and parameter t (variance).
%                  This routine does not employ the iterative solver. Instead,
%                  it takes advantage of the roots solve of matlab.
%                  The algorithm is much faster and stable for high order 
%                  Hermite polynomials.
%
% Syntax :    z = HermiteZeros_t( degree, t )
%
% Input  :    degree - degree of Hermite polynomial (number of zeros)
%             t > 0 is the parameter.
%
% Output :    z - zeros as array (degree x 1)
%
% By Dongbin Xiu   12/12/2002
%

z = roots(HermiteCoef_t(degree,t));
