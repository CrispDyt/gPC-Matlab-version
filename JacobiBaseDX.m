function dbase = JacobiBaseDX(x, P)
%
% JacobBaseDX.m - Compute the derivatives Jacobi polynomial basis functions
%                 in space.
%
% Syntax:   based = JacobBaseDX(x, mode);
%
% Input :   x = x-coordinate in row vector form -1<x<1  
%           P = the highest number of mode 0<= p <=P
%
% Output:   dbase = bases function derivatives stored at location x stored as
% (P+1) x length(x).
%
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% Written by Dongbin Xiu   3/26/2004. 
%

np = length(x);

dbase = zeros(P+1,np);

dbase(1,:) = -0.5;
dbase(P+1,:) = 0.5;
for p=1:P-1
  dbase(p+1,:) = -0.5*(1+x)/2.*JacobiF(x, p-1, 1, 1) + ...
                  0.5*(1-x)/2.*JacobiF(x, p-1, 1, 1) + ...
                (1-x)/2.*(1+x)/2.*JacobiD(x, p-1, 1, 1);
end
