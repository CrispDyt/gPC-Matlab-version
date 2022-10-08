function z = CharlierZeros(n, alpha)
%
% CharlierZeros.m - Evaluates the zeros of Charlier polynomial
%                  with degree n and parameter alpha.
%                  This routine does not employ the iterative solver. Instead,
%                  it takes advantage of the roots solve of matlab.
%                  the limit of the algorithm is reached.
%
% Syntax :    z = CharlierZeros( n, alpha )
%
% Input  :    n - degree of Charlier polynomial (number of zeros)
%
% Output :    z - zeros as column vector (n x 1) in ascending order.
%
% By Dongbin Xiu   5/02/2003
%

z = sort(roots(CharlierCoef(n, alpha)));
