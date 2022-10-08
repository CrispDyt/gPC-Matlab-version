function base = JacobiBaseX(x, P)
%
% JacobF.m - Compute the Jacobi polynomial basis function
%            in space.
%
% Syntax:   base = jacobf(x, mode, alpha, beta);
%
% Input :   x = x-coordinate in row vector form -1<x<1  
%           P = the highest number of mode 0<= p <=P
%
% Output:   base = bases function stored at location x stored as
% (P+1) x length(x).
%
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% Written by Dongbin Xiu   10/21/2003. Tested by evaluating the
% mass matrix structure (see also JacobiMass.m).
%

np = length(x);

base = zeros(P+1,np);

base(1,:) = (1-x)/2;
base(P+1,:) = (1+x)/2;
for p=1:P-1
  base(p+1,:) = (1-x)/2.*(1+x)/2.*JacobiF(x, p-1, 1, 1);
end
