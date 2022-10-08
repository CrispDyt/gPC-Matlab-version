function coef = HermiteCoef( n )
%
% HermiteCoef.m - returns the coefficients of Hermite polynomials in monic form
%                 H^n(x) = x^0 + a_1 x + ... + a_n x^n, where a_n = 1;
%                 This is the format in matlab for polynomials.
%
% Syntax :   coef = HermiteCoef(n)
%
% Input :    n - order of the Hermite polynomial
% Output:    coef - (1 x (n+1)) vector containing a_i, the coefficients of 
%            x^(n-i), where i=0,...n. Note half of the coefficients are zero.
%
% Results were validated against analytical formula for degree (n) up to 50,
% where the larger coefficients are of the order 10^33.
%
% by Dongbin Xiu  5/08/2002. Modified the output format 12/09/12.
%

coef = zeros(1,n+1);

for r = 0:round((n-0.1)/2)
  coef(n-2*r+1) = (-1)^r/2^r * nchoosek(n,r) * prod((n-2*r+1):(n-r));
end

coef=fliplr(coef);
