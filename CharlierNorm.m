function d = CharlierNorm(n, alpha)
%
% CharlierNorm.m - Evaluates the square norms of Charlier polynomial
%                  with degree n and parameter alpha.
%                  d^2 = alpha^n * n!
%
% Syntax :    d = CharlierNorm( n, alpha )
%
% Input  :    n - degree of Charlier polynomial.
%
% Output :    d - the square norm of n-degree Charlier polynomial
%
% By Dongbin Xiu   5/05/2003
%

d = alpha^n * factorial(n);
