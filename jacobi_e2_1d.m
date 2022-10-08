function e = jacobi_e2_1d(order, alpha, beta)
%
% jacobi_e2_1d.m - Evaluate the inner product of 1d Jacobi-chaos doubles
%
% Syntax     e = jacobi_e2_1d(order, alpha, beta)
%
% Input:     order = order of Jacobi-chaos
%            alpha, beta = parameters of Jacobi-chaos (alpha, beta>-1)
% Output:    e = 1 x (p+1) row vector containing the result.
% 
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% By Dongbin Xiu   04/13/2003
%

tolerance=1e-10;
e=zeros(order+1,1);

np = ceil((2*order+1)/2);

[z,w] = JacobiZW(np, alpha, beta);  

factor = 2^(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)/gamma(alpha+beta+2);

for i=1:order+1
   e(i) = sum(JacobiF(z,i-1,alpha,beta).^2.*w)/factor;
end



