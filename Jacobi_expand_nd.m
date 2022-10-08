function f = Jacobi_expand_nd(x, coef, alpha, beta, norder)
%
% Jacobi_expand_nd.m - Compute all multi-dimensional Jacobi polynomials
%           with a given order at given coordinates in arbitrary dimension.
%
% Syntax:   f = Jacobi_expand_nd(x, ndim, norder, alpha, beta);
%
% Input :   x = (npt,ndim) array containing 'npt' number ndim-coordinates
%           ndim = dimension of space
%           order = highest order of polynomials returned.
%           alpha, beta = parameters of Jacobi polynomial (alpha,beta>-1)
% Output:   f = (npt, 1) array. Polynomials evaluated at x.
%
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% by Dongbin Xiu   10/09/2007
%

npt = length(x(:,1));
ndim = length(x(1,:));

%pmatrix = chaos_sequence(ndim,norder);
nterm = min(nchoosek(ndim+norder,ndim), length(coef));

f=zeros(npt,1);

for m=1:npt
  poly = JacobiF_nd(x(m,:), ndim, norder, alpha, beta);
  for n=1:nterm
    f(m) = f(m) + coef(n)*poly(n);
  end 
end
