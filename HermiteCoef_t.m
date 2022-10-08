function coef = HermiteCoef_t( n, t )
%
% HermiteCoef_t.m - returns the coefficients of Hermite polynomials as
%                 H^n(x,t) = x^n + a_(n-1) x^(n-1) + ... + a_0, where a_n = 1;
%                 The coefficients are stored in matlab format.
%
% Syntax :   coef = HermiteCoef_t(n,t)
%
% Input :    n - order of the Hermite polynomial
% Output:    coef - (1 x (n+1)) vector containing a_i, the coefficients of x^i,
%            where i=0,...n. Note half of the coefficients are zero.
%
% Results were validated indirectly for the eigenvalue problem from Hermite
% chaos expansion, where the result matches for degree up to 45.
%
% by Dongbin Xiu  11/21/2002. Change of output format made 12/12/02.
%

coef = HermiteCoef(n);

cnt=0;
for r = 1:2:(n+1)
  coef(r) = coef(r)*(t^cnt);
  cnt = cnt + 1;
end

