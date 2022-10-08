function e = Laguerre_e2_nd( ndim, p, alpha)
%
% Laguerre_e2_nd.m - Evaluate the inner product of n-dimensional
%                  Laguerre-chaos, i.e. int L^2 dx
%
% Syntax     e = Laguerre_e2_nd( ndim, p, alpha)
%
% Input:     ndim = dimensionality of the Jacobi-chaos
%            p    = order of the Jacobi-chaos
%            alpha = parameters of Jacobi-chaos (alpha>-1)
% Output:    e = 1xP (P=(ndim+p)!/(ndim!p!)) array containing the result.
% 
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% By Dongbin Xiu   05/22/2007
%

if ndim <= 0
   fprintf('Invalid dimensionality! (need n>0)\n');
end

if ndim == 1
   e = Laguerre_e2_1d(p, alpha);
else
   e1 = Laguerre_e2_1d(p, alpha);
   poly = chaos_sequence(ndim, p);
   [P,ntmp] = size(poly);
   e = ones(P,1);
   for m=1:P
       for n=1:ndim
           e(m) = e(m) * e1(poly(m,n)+1);
       end
   end
end

