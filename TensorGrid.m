function z = TensorGrid(z1, ndim)
%
% The routine takes a one-dimensional vector and construct ndim-dimensional
% tensor grid.
%
% z = TensorGrid(z1, ndim)
% 
% z1 -- vector of length np, row or column vector.
% ndim -- dimensionality
% z -- ndim-dimensional tensor grid of size (ntot, ndim),
%      where ntot = np^ndim.
%
% by Dongbin Xiu, 11/17/2005.
%

np = length(z1);
if ndim == 1
    z = zeros(np,ndim);
    z(:) = z1(:);
elseif ndim == 2
    z = [z1(1) z1(1)];
    for i=1:np
        for j=1:np
            if  ((i~=1) | (j~=1))
                z = [z; [z1(i) z1(j)]];
            end
        end
    end
else
    z = TensorGrid(z1,ndim-1);
end
