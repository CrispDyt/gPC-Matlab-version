function e = jacobi_e2_nd( ndim, p, alpha, beta)
%
% jacobi_e2_nd.m - Evaluate the inner product of n-dimensional
%                  Jacobi-chaos, i.e. int J^2 dx
%
% Syntax     e = jacobi_e2_nd( ndim, p, alpha, beta)
%
% Input:     ndim = dimensionality of the Jacobi-chaos
%            p    = order of the Jacobi-chaos
%            alpha, beta = parameters of Jacobi-chaos (alpha, beta>-1)
% Output:    e = 1xP (P=(ndim+p)!/(ndim!p!)) array containing the result.
% 
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% By Dongbin Xiu   11/09/2005
%

if ndim <= 0
   fprintf('Invalid dimensionality! (need n>0)\n');
end

if ndim == 1
   e = jacobi_e2_1d(p, alpha, beta);
else
   e1 = jacobi_e2_1d(p, alpha, beta);
   poly = chaos_sequence(ndim, p);
   [P,ntmp] = size(poly);
   e = ones(P,1);
   for m=1:P
       for n=1:ndim
           e(m) = e(m) * e1(poly(m,n)+1);
       end
   end
end

