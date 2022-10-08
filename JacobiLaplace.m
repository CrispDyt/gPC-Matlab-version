function M = JacobiLaplace(P)
%
% JacobiLaplace.m - Compute the stiffness matrix of the Laplace operator
%                   in 1D and its inverse from
%                   Jacobi polynomial basis function in space.
%
% Syntax:   M = JacobiLaplace(P)
%
% Input :   P = the highest number of mode 0<= p <=P
%
% Output:   M are (P+1)x(P+1) matrices, mass matrix and its inverse.
%
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% Written by Dongbin Xiu   3/26/2004. 
% Stiff matrix structure examined. (See book by GK and Spencer, pp 49.)

M = zeros(P+1,P+1);
np = ceil((2*P+1)/2);
alpha=0;   beta=0;
factor = 2^(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)/gamma(alpha+ ...
                                                  beta+2);

[z,w] = JacobiZW(np,alpha,beta);

bp = JacobiBaseDX(z',P);
for i=1:P+1
  for j=1:P+1
    M(i,j) = sum(bp(i,:).*bp(j,:).*w')/factor;
    if abs(M(i,j)) < 1e-12
      M(i,j) = 0;
    end
  end
end

