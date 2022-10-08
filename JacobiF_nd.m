function poly = JacobiF_nd(x, ndim, norder, alpha, beta)
%
% JacobiF_nd.m - Compute all multi-dimensional Jacobi polynomials
%           up to a given order at given coordinates.
%
% Syntax:   poly = JacobiF_nd(x, ndim, norder, alpha, beta);
%
% Input :   x = (1,ndim) array containing a ndim-coordinates
%           ndim = dimension of space
%           norder = highest order of polynomials returned.
%           alpha, beta = parameters of Jacobi polynomial (alpha,beta>-1)
% Output:   poly = (1, nterm) array. All nterm polynomials evaluated at x.
%
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% by Dongbin Xiu   11/08/2005
%

ntmp = length(x);

if ntmp ~= ndim
fprintf('Incompatible dimension in JacobiF_nd! Quit.\n');
end

pmatrix = chaos_sequence(ndim,norder);
ptmp = zeros(norder+1,ndim);
for n=0:norder
    ptmp(n+1,:) = JacobiF(x,n,alpha,beta);
end

[nterm,ntmp]=size(pmatrix);
poly = ones(1, nterm);

for m=2:nterm
    for n=1:ndim
        poly(m) = poly(m) * ptmp(pmatrix(m,n)+1,n);
    end
end