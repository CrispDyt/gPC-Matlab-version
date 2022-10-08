function e2 = hermite_ij_nd( ndim, norder )
%
% The function returns the double inner products of Hermite-chaos basis
% of order norder and dimension ndim.
%
% e2 - array (nterm x 1) containing the results, i.e., the diagonal terms
% of the diagonal matrix. Here nterm=nchoosek(ndim_norder, ndim) is the
% total number of bases.
%
% Written by Dongbin Xiu 2/24/2004.
%
% To evaluate the inner product of an arbitrary set of basis, use
% hermite_h2_nd.
%

ip = zeros(1,ndim);

index = index_step( ip, norder );

[nterm, ntmp] = size( index );
e2 = zeros( nterm, 1 );

for i = 1:nterm
   e2(i) = hermite_h2( index(i,:), index(i,:));
end
