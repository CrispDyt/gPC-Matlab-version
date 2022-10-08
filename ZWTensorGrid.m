function [z,w] = ZWTensorGrid(z1, w1, ndim)
%
% The routine takes a one-dimensional vector and construct ndim-dimensional
% tensor grid.
%
% z = TensorGrid(z1, w1, ndim)
% 
% z1 -- vector of length np, row or column vector. 1d grid.
% w1 -- 1d weights.
% ndim -- dimensionality
% z -- ndim-dimensional tensor grid of size (ntot, ndim),
%      where ntot = np^ndim.
%
% by Dongbin Xiu, 5/12/2006.
%

np = length(z1);
if ndim == 1
    z = zeros(np,ndim);  w=zeros(np,1);
    z(:) = z1(:);  w(:)=w1(:);
elseif ndim == 2
    z = [z1(1) z1(1)];   w=w1(1)*w1(1);
    for i=1:np
        for j=1:np
            if  ((i~=1) | (j~=1))
                z = [z; [z1(i) z1(j)]];
                w = [w; w1(i)*w1(j)];
            end
        end
    end
else
    [z,w] = ZWTensorGrid(z1,w1,ndim-1);
end