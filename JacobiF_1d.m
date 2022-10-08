function poly_out = JacobiF_1d(x, n, alpha, beta)
%
% JacobiF.m - Compute the Jacobi polynomial P^n(x; alpha, beta)
%
% Syntax:   poly_out = JacobiF_1d(x, n, alpha, beta);
%
% Input :   x = x-coordinate of one point
%           n = order of polynomial
%           alpha, beta = parameters of Jacobi polynomial (alpha,beta>-1)
% Output:   poly_out = polynomial values at location x of orders 0 to n.
%
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% by Dongbin Xiu   06/272011
%

one = 1;
two = 2;
np=1;

poly_out = ones(1,n+1);


if n == 0
poly = ones(np);
poly_out(1) = 1;

elseif  n == 1
poly = 0.5*(alpha - beta + (alpha + beta + two)*x);
poly_out(n+1) = poly;

else
apb = alpha + beta;
    
polyn2 = ones(np);
polyn1 = 0.5*(alpha - beta + (alpha + beta + two)*x);

    
for k = 2:n
	if k == 2
		poly_out(k) = polyn1;
	end
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

	poly_out(k+1) = poly;
end




end  
