function [M, IM] = JacobiMass(P)
%
% JacobiMass.m - Compute the mass matrix and its inverse from
%                Jacobi polynomial basis function in space.
%
% Syntax:   [M,IM] = JacobiMass(P)
%
% Input :   P = the highest number of mode 0<= p <=P
%
% Output:   M, IM are (P+1)x(P+1) matrices, mass matrix and its inverse.
%
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% Written by Dongbin Xiu   10/21/2003. 
% Mass matrix structure examined. (See book by GK and Spencer.)

M = zeros(P+1,P+1);  IM = M;
np = ceil((2*P+1)/2);
alpha=0;   beta=0;
factor = 2^(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)/gamma(alpha+ ...
                                                  beta+2);

[z,w] = JacobiZW(np,alpha,beta);

bp = JacobiBaseX(z',P);
for i=1:P+1
  for j=1:P+1
    M(i,j) = sum(bp(i,:).*bp(j,:).*w')/factor;
    if abs(M(i,j)) < 1e-12
      M(i,j) = 0;
    end
  end
end

IM=inv(M);