function coef = LaguerreCoef( n, alpha )
%
% LaguerreCoef.m - returns the coefficients of Laguerre polynomials in 
%                  monic form,
%                  L^n(x,alpha) = a_1 x^n + a_2 x^n-1 + ... + a_n x + a_{n+1};
%                  This is the format in matlab for polynomials.
%
% Syntax :   coef = LaguerreCoef(n, alpha)
%
% Input :    n - order of the Laguerre polynomial L^n(x,alpha)
% Output:    coef - (1 x (n+1)) vector containing a_i, the coefficients of 
%            x^(n-i), where i=0,...n. 
%
% by Dongbin Xiu  12/12/2002. 
%

coef = zeros(1,n+1);

for k = 0:n
  coef(k+1) = pochhammer(-n,k)*pochhammer(alpha+k+1,n-k)/factorial(k);
end

coef=fliplr(coef/factorial(n));
