function polyd = jacobd(x, n, alpha, beta)
%
% jacobd.m - Evaluates the derivatives of jacobi polynomials P^n(x).
%  
% Syntax:   polyd = jacobd(x, n, alpha, beta);
%
% Input :   x = points where the derivatives are to be evaluated
%           (n, alpha, beta) = parameters of the Jacobi polynomial
%           P^n(x; alpha, beta).
% Output:   polyd = values of derivatives at x, stored in same form as x.
%
%  To do this we have used the relation 
%
%  d   alpha,beta   1                  alpha+1,beta+1
%  -- P (z)       = -(alpha+beta+n+1) P  (z)
%  dz  n            2                  n-1
%
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% Exported from Nektar library by Dongbin Xiu. 01/23/2002.
%

np=size(x);

if n == 0
polyd = zeros(np);

else
polyd = 0.5*(alpha+beta+n+1)*jacobf(x,n-1,alpha+1,beta+1);

end


