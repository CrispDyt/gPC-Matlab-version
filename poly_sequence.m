function [q, kmat, fac_mat] = poly_sequence(ndim, p)
%
% poly_sequence.m -- returns the order matrix for multi-dimensional
%                     polynomials.
%
% Syntax:  q = poly_sequence(ndim,p)
%          [q, kmat] = poly_sequence(ndim,p)
%
% Input: ndim - dimensionality
%        p    - highest order
%
% Output: q - matrix of size (nterm, ndim)
%         where nterm = nchoosek(ndim+p,ndim)
%         kmat - matrix of size (p,2), where each row 'k' stores
%                the location [k_low, k_upp] indicating the index range
%                covering the k-th order terms.
%

if ndim <= 0
   fprintf('Invalid dimensionality! (need n>0)\n');
else
   %q=index_step(zeros(1,ndim),p); 
   kmat = zeros(p,2);
   
   % zero degree
   q = zeros(1,ndim);   
   kmat(1,1) = 1;  kmat(1,2) = 1; 

   qtmp = q;
   for i=1:p
      q1 = index_step1(qtmp);
      [nlength, nd] = size(q1);
      q = cat(1, q, q1);
      kmat(i+1,1) = kmat(i,2) + 1;   
      kmat(i+1,2) = kmat(i,2) + nlength;
      qtmp = q1;
   end

   if nargout == 3    % calculate the multiindex factorials
     [nterm, ndim] = size(q);
     fac_mat = ones(nterm,1);
     for i=2:nterm
       for j=1:ndim
         fac_mat(i) = fac_mat(i) * factorial(q(i,j));
       end
     end
   end
       
end
