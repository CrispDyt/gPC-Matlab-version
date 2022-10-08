function q = chaos_sequence(ndim, p)
%
% chaos_sequence.m -- returns the order matrix for multi-dimensional
%                     polynomials.
%
% Syntax:  q = chaos_sequence(ndim,p)
%
% Input: ndim - dimensionality
%        p    - highest order
%
% Output: q - matrix of size (nterm, ndim)
%         where nterm = nchoosek(ndim+p,ndim)
%

if ndim <= 0
   fprintf('Invalid dimensionality! (need n>0)\n');
else
   q=index_step(zeros(1,ndim),p);
end
