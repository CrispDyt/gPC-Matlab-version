function poly = HermiteD_nd(x, ndim, norder, nder)
%
% HermiteD_nd.m - Compute 1st derivative of multi-dimensional Hermite polynomials
%           up to a given order at a given coordinate.
%
% Syntax:   poly = HermiteD_nd(x, ndim, norder, nder);
%
% Input :   x = (1,ndim) array containing a ndim-coordinates
%           ndim = dimension of space
%           norder = highest order of polynomials returned.
%           nder = the dimensional where derivative is taken.
%
% Output:   poly = (1, nterm) array. All nterm polynomials evaluated at x.
%
% by Dongbin Xiu   11/17/2005
%

ntmp = length(x);

if ntmp ~= ndim
    fprintf('Incompatible dimension in JacobiF_nd! Quit.\n');
end

pmatrix = chaos_sequence(ndim,norder);
ptmp = zeros(norder+1,ndim);
for n=0:norder
    ptmp(n+1,:) = HermiteF(x,n);
    pder(n+1,:) = HermiteD(x,n);
end

[nterm,ntmp]=size(pmatrix);
poly = ones(1, nterm);
poly(1) = 0;

for m=2:nterm
    for n=1:ndim
        if n == nder
            poly(m) = poly(m) * pder(pmatrix(m,n)+1,n);
        else
            poly(m) = poly(m) * ptmp(pmatrix(m,n)+1,n);
        end
    end
end