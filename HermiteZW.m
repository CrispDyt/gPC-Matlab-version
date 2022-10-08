function [z,w] = HermiteZW( degree )
%
% HermiteZW.m - Evaluates the zeros and weights of Hermite polynomial
%               with degree 'degree'.
%
% Syntax :    [z,w] = HermiteZW( degree )
%
% Input  :    degree - degree of Hermite polynomial (number of zeros)
%
% Output :    [z, w] - zeros and weights in array (degree x 1)
%
% By Dongbin Xiu   5/03/2002
%

nfac = factorial(degree);
twonp1 = 2^(degree+1);

%z = HermiteZeros(degree);  
z = HermiteZeros_iter(degree);

w = sqrt(pi)*nfac./(HermiteD(z,degree).^2);
