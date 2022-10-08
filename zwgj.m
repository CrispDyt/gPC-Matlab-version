function [z,w] = zwgj(np, alpha, beta)
%
% zwgj.m  -  Evaluates the Gauss-Jacobi quadrature points and the corresponding
%            weights. The quadrature points are the roots of P_np(alpha,beta).
%
% Syntax:  [z,w] = zwgj(np, alpha, beta)
%
% Input :  np = number of quadrature points
%          alpha, beta = parameters define the Jacobi polynomial, alpha,beta>-1
%
% Output:  [z,w] = quadrature points and weights stored in array with length np
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% Exported from Nektar library by Dongbin Xiu   01/24/2002
%

apb = alpha + beta;

z = jacobz(np, alpha, beta);
w = jacobd(z, np, alpha, beta);

fac = 2^(apb+1)*gamma(alpha+np+1)*gamma(beta+np+1);
fac = fac/(gamma(np+1)*gamma(apb+np+1));
  
w = fac./((w.^2).*(1-z.*z));

