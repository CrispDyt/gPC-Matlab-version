function [z,w] = CharlierZW( n, alpha )
%
% CharlierZW.m - Evaluates the zeros and weights of Charlier polynomial
%               with degree n and parameter alpha.
%
% Syntax :    [z,w] = CharlierZW( n, alpha )
%
% Input  :    n - degree of Charlier polynomial,
%             alpha is the parameter
%
% Output :    [z, w] - zeros and weights in column vector (n x 2)
%
% By Dongbin Xiu   5/05/2003
%

z = CharlierZeros(n, alpha);
d = CharlierNorm(n-1,alpha);

f  = CharlierCoef(n-1, alpha);
df = polyder(CharlierCoef(n, alpha));

w = d./(polyval(df, z).*polyval(f, z));
