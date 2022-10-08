function z = JacobiZeros_direct(degree, alpha, beta)
%
% JacobiZeros.m - Evaluates the zeros of Jacobi polynomials
%                  with degree 'degree'.
%                  This routine does not employ the iterative solver. Instead,
%                  it takes advantage of the roots solve of matlab.
%                  The algorithm is slower and less stable for high order 
%                  Jacobi polynomials, in contrast to the similar Hermite
%                  algorithm.
%
% Syntax :    z = JacobiZeros_direct( degree, alpha, beta )
%
% Input  :    degree - degree of Jacobi polynomial (number of zeros) with
%                      parameters alpha, beta > -1.
%
% Output :    z - zeros as array (degree x 1)
%
% By Dongbin Xiu   12/09/2002
%

z = roots(JacobiCoef(degree, alpha, beta));
