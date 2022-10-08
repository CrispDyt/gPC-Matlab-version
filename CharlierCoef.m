function coef = CharlierCoef(n, alpha)
%
% CharlierF.m - Compute the coefficients of Charlier polynomial C^n(x;alpha)
%               in monic form,
%               C^n(x,alpha) = a_1 x^n + a_2 x^n-1 + ... + a_n x + a_{n+1};
%               This is the format in matlab for polynomials.
%               Computed Charlier polynomials have unit leading coefficients.
%
% Syntax:   poly = CharlierF(x, n, degree);
%
% Input :   x = x-coordinate in matrix form (integers)
%           n = order of polynomial,
%           alpha > 0 is the parameter.
%
% Output:   poly = polynomial values at location x stored in same form as x.
%
% This routine uses a formula recursively, but not the recurrence formula.
%
% by Dongbin Xiu   5/05/2003.
%
coef = zeros(1,n+1);

for k=0:n
xkcoef=[fliplr(xchoosek_coef(k)) zeros(1,n-k)];
coef=coef+xkcoef*nchoosek(n,k)*factorial(k)*(-alpha)^(n-k);
end

coef=fliplr(coef);

function xkcoef=xchoosek_coef(k)

xkcoef = zeros(1,k+1);
switch k
   case 0
      xkcoef=1;
   otherwise
      xklast = xchoosek_coef(k-1);
      xkcoef = ([xklast 0] - (k-1)*[0 xklast])/k;
end

