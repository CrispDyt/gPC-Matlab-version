function poly = JacobiF(x, n, alpha, beta)
%
% JacobiF.m - Compute the Jacobi polynomial P^n(x; alpha, beta)
%
% Syntax:   poly = JacobiF(x, n, alpha, beta);
%
% Input :   x = x-coordinate in matrix form
%           n = order of polynomial
%           alpha, beta = parameters of Jacobi polynomial (alpha,beta>-1)
% Output:   poly = polynomial values at location x stored in same form as x.
%
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% Exported from Nektar library by Dongbin Xiu   01/23/2002
%

np = size(x);
one = 1;
two = 2;

if n == 0
poly = ones(np);

elseif  n == 1
poly = 0.5*(alpha - beta + (alpha + beta + two)*x);

else
apb = alpha + beta;
    
polyn2 = ones(np);
polyn1 = 0.5*(alpha - beta + (alpha + beta + two)*x);
    
for k = 2:n
a1 =  two*k*(k + apb)*(two*k + apb - two);
a2 = (two*k + apb - one)*(alpha^2 - beta^2);
a3 = (two*k + apb - two)*(two*k + apb - one)*(two*k + apb);
a4 =  two*(k + alpha - one)*(k + beta - one)*(two*k + apb);
      
a2 = a2/a1;
a3 = a3/a1;
a4 = a4/a1;
	
poly = (a2 + a3*x).*polyn1 - a4*polyn2;
polyn2 = polyn1;
polyn1 = poly;
end

end  
