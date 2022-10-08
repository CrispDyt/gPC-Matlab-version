function z = HermiteZeros(degree)
%
% HermiteZeros.m - Evaluates the zeros of Hermite polynomial
%                  with degree 'degree'.
%                  This routine does not employ the iterative solver. Instead,
%                  it takes advantage of the roots solve of matlab.
%                  The algorithm is much faster and stable for high order 
%                  Hermite polynomials. Accuracy is ensured for order lower 
%                  than 50; however orders less than 60 are still considered
%                  good despite the loss of accuracy from binomial coefficients.
%                  For orders greater than 60, complex roots appear indicating
%                  the limit of the algorithm is reached.
%
% Syntax :    z = HermiteZeros( degree )
%
% Input  :    degree - degree of Hermite polynomial (number of zeros)
%
% Output :    z - zeros as array (degree x 1)
%
% By Dongbin Xiu   12/09/2002
%

z = roots(HermiteCoef(degree));
