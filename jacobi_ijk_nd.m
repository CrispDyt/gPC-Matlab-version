function e = jacobi_ijk_nd( ndim, p, alpha, beta)
%
% jacobi_ijk_nd.m - Evaluate the inner product of n-dimensional
%                  Jacobi-chaos triplets
%
% Syntax     e = jacobi_ijk_nd( ndim, p, alpha, beta)
%
% Input:     ndim = dimensionality of the Jacobi-chaos
%            p    = order of the Jacobi-chaos
%            alpha, beta = parameters of Jacobi-chaos (alpha, beta>-1)
% Output:    e = PxPxP (P=(ndim+p)!/(ndim!p!)) matrix containing the result.
% 
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% By Dongbin Xiu   04/17/2002
%

if ndim <= 0
   fprintf('Invalid dimensionality! (need n>0)\n');
end

if ndim == 1
   e = jacobi_ijk_1d(p, alpha, beta);
else
   P = nchoosek(ndim+p, p);
   e = ones(P,P,P);
   e1 = jacobi_ijk_1d(p, alpha, beta);
   poly = chaos_sequence(ndim, p);
   for i=1:P
      for j=i:P
         for k=j:P
            for n=1:ndim
               e(i,j,k) = e(i,j,k) * e1(poly(i,n)+1, poly(j,n)+1, poly(k,n)+1);
            end

            e(i,k,j) = e(i,j,k);
            e(j,i,k) = e(i,j,k);
            e(j,k,i) = e(i,j,k);
            e(k,i,j) = e(i,j,k);
            e(k,j,i) = e(i,j,k);
         end
      end
   end
end

