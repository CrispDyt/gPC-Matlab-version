function poly = LaguerreF_nd(x, ndim, norder, alpha)
%
% LaguerreF_nd.m - Compute multi-dimensional Hermite polynomials
%           up to a given order at a given coordinate.
%
% Syntax:   poly = LaguerreF_nd(x, ndim, norder, alpha);
%
% Input :   x = (1,ndim) array containing a ndim-coordinates
%           ndim = dimension of space
%           norder = highest order of polynomials returned
%           alpha  = real parameter > -1.
%
% Output:   poly = (1, nterm) array. All nterm polynomials evaluated at x.
%
% by Dongbin Xiu   5/07/2007
%

ntmp = length(x);

if ntmp ~= ndim
    fprintf('Incompatible dimension in JacobiF_nd! Quit.\n');
end

pmatrix = chaos_sequence(ndim,norder);
ptmp = zeros(norder+1,ndim);
for n=0:norder
    ptmp(n+1,:) = LaguerreF(x,n,alpha);
end

[nterm,ntmp]=size(pmatrix);
poly = ones(1, nterm);

for m=2:nterm
    for n=1:ndim
        poly(m) = poly(m) * ptmp(pmatrix(m,n)+1,n);
    end
end