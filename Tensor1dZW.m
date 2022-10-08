function [z,w] = Tensor1dZW(z1,w1, ndim)
%
% The routine takes a one-dimensional quadrature rule and its n-dim
% counterpart.
%
% z = Tensor1dZW(z1, w1, ndim)
% 
% z1, w1 -- vector of length Np, row or column vector for 1d points +
% weights.
% ndim -- dimensionality
% z -- ndim-dimensional tensor grid of size (ntot, ndim),
%      where ntot = Np^ndim.
% If ndim==1, then [z,w] returns the same but as column vectors.
%
% by Dongbin Xiu, 11/18/2013.
%

[nc, nr] = size(z1);

if nc == 1
    z1 = z1';  w1 = w1';
end
    
if ndim == 1
    z = z1; w = w1;
else
    z = z1; w=w1;
    for d=2:ndim
        [z,w] = TensorZW2(z,w,z1,w1);
    end
end
