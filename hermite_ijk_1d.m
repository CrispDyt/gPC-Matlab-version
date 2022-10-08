function e = hermite_ijk_1d_direct(p)
%
% hermite_ijk_1d.m - Evaluate the inner product of 1d Hermite-chaos triplets
%     via the exact formula, which is faster than quadrature at high orders.
%
% Syntax     e = hermite_ijk_1d_direct(p)
%
% Input:     p = order of Hermite-chaos
%
% Output:    e = (p+1)x(p+1)x(p+1) matrix containing the result as integers.
% 
% By Dongbin Xiu   11/19/2002
%

e = zeros(p+1, p+1, p+1);
  
for i=1:p+1
   for j=i:p+1
      for k=j:p+1
         s = i + j + k - 3;
         if mod(s,2) == 0
            s = s / 2;
            if ((s >= (i-1)) & (s>=(j-1)) & (s>=(k-1)))
	       e(i,j,k) = factorial(i-1)*factorial(j-1)*factorial(k-1)/ ...
                       (factorial(s-i+1)*factorial(s-j+1)*factorial(s-k+1));
 
	       e(i,k,j) = e(i,j,k);
	       e(j,i,k) = e(i,j,k);
	       e(j,k,i) = e(i,j,k);
	       e(k,i,j) = e(i,j,k);
	       e(k,j,i) = e(i,j,k);
            end
         end
      end
   end
end
