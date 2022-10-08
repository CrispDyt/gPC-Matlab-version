function e2 = hermite_h2_nd( index )
% 
% This function evaluates the inner products of Hermite-chaos basis
% <H_i H_j> = a^2 \delta_ij. The exact formula is employed.
%
% index - array of size (nterm x ndim), where nterm is the total number
% of bases involved (arbitrary set of basis) with dimension ndim.
% Each row represents a base where each entry of the row vector is the
% order of Hermite polynomial in that dimension.
%
% e2 - (nterm x 1) array containing the results. Note formally the result
% should be a diagonal matrix.
%
% Written by Dongbin Xiu 2/24/2004.
%

[nterm, ndim] = size( index );
e2 = zeros( nterm, 1 );

for i = 1:nterm
   e2(i) = hermite_h2( index(i,:), index(i,:));
end
