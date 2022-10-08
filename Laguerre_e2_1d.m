function e = Laguerre_e2_1d(order, alpha)
%
% Laguerre_e2_1d.m - Evaluate the inner product of 1d Laguerre-chaos doubles
%
% Syntax     e = Laguerre_e2_1d(order, alpha)
%
% Input:     order = order of Laguerre-chaos
%            alpha = parameters of Laguerre-chaos (alpha>-1)
% Output:    e = 1 x (p+1) row vector containing the result.
% 
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% By Dongbin Xiu   05/22/2007
%

e=zeros(order+1,1);


factor = gamma(alpha+1);   % factor to normaliz constant integral

for n=0:order
   e(n+1) = (gamma(n+alpha+1)/factorial(n))/factor;
end



