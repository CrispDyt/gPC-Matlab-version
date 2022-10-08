function poly = JacobiF_scale(x, n, alpha, beta, s)
%
% JacobF.m - Compute the scaled Jacobi polynomial P^n(x; alpha, beta, s)
%
% Syntax:   poly = jacobf(x, n, alpha, beta, s);
%
% Input :   x = x-coordinate in matrix form
%           n = order of polynomial
%           alpha, beta = parameters of Jacobi polynomial (alpha,beta>-1)
%           s = the scaling factor, s>0;
% Output:   poly = polynomial values at location x stored in same form as x.
%
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% Exported from Nektar library by Dongbin Xiu   09/24/2004
%

if s==1
  poly = JacobiF(x,n,alpha,beta);
else
  np = size(x);
  one = 1;
  two = 2;

  x = s*x;

  if n == 0
     poly = ones(np);
  elseif  n == 1
     poly = 0.5*(alpha - beta + (alpha + beta + two)*x);
  else
   apb = alpha + beta;

   k=n;
   a1 =  two*k*(k + apb)*(two*k + apb - two);
   a2 = (two*k + apb - one)*(alpha^2 - beta^2);
   a3 = (two*k + apb - two)*(two*k + apb - one)*(two*k + apb);
   a4 =  two*(k + alpha - one)*(k + beta - one)*(two*k + apb);
      
   a2 = a2/a1;
   a3 = a3/a1;
   a4 = a4/a1;

   poly = (a2*s + a3*x).*JacobiF_scale(x,n-1,alpha,beta,s) - ...
          s^2*a4*JacobiF_scale(x,n-2,alpha,beta,s);

end  
end
