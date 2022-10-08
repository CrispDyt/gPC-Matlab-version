function q = hyperbolic_sequence(ndim, p)
%
% hyperbolic_sequence.m -- returns the order matrix for multi-dimensional
%                     polynomials.
%
% Syntax:  q = hyperbolic_sequence(ndim,p)
%
% Input: ndim - dimensionality
%        p    - highest order
%
% Output: q - matrix of size (nterm, ndim)
%         where nterm is the size of all indices satisfying hyperbolic
%         cross
%

if ndim <= 0
   fprintf('Invalid dimensionality! (need n>0)\n');
else
   q=index_step_hc(zeros(1,ndim),p);
end
