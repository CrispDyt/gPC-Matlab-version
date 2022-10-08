function z = jacobz(n, alpha, beta)
%
% jacobz.m  -  Evaluates the zeros of n-th order Jacobi polynomial
%
% Syntax:  z = jacobz(n, alpha, beta)
% 
% Input :  n = order of Jacobi polynomial P^n(x; alpha, beta)
%          alpha, beta = parameters define the polynomial, alpha,beta>-1
%
% Output:  z = zeros of the polynomial, stored in array with length of n.
%
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% Exported from Nektar library by Dongbin Xiu   01/24/2002
%

STOP=100;    % maximum number of iteration
EPS =1e-12;  % error tolerance of iteration

dth = pi/(2*n);
rlast=0;

z = zeros(n,1);
  
for k = 1:n
r = -cos((2*k - 1) * dth);
if k>1
r = (r + rlast)/2;
end
    
for j = 1:STOP
poly = jacobf(r, n, alpha, beta);
pder = jacobd(r, n, alpha, beta);
    
sum=0;  
for i = 1:(k-1) 
sum = sum+ 1/(r - z(i));
end

delr = -poly(1) / (pder(1) - sum * poly(1));
r    = r + delr;

if abs(delr) < EPS  
break;
end

end

z(k)  = r;
rlast = r;

end

z=sort(z);
