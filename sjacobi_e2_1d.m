function e = sjacobi_e2_1d(order, alpha, beta, s)
%
% jacobi_e2_1d.m - Evaluate the inner product of 1d scaled Jacobi-chaos doubles
%
% Syntax     e = sjacobi_e2_1d(order, alpha, beta,s)
%
% Input:     order = order of Jacobi-chaos
%            alpha, beta = parameters of Jacobi-chaos (alpha, beta>-1)
% Output:    e = 1 x (p+1) row vector containing the result.
% 
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% By Dongbin Xiu   09/24/2004
%

tolerance=1e-10;
e=zeros(1,order+1);

np = ceil((2*order+1)/2);

[z,w] = zwgj(np, alpha, beta);  

factor = 2^(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)/gamma(alpha+beta+2);

for i=1:order+1
   e(i) = sum(JacobiF_scale(z,i-1,alpha,beta,s).^2.*w)/factor;
end



