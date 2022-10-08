function e = hermite_ijk_1d(p)
%
% hermite_ijk_1d.m - Evaluate the inner product of 1d Hermite-chaos triplets
%
% Syntax     e = hermite_ijk_1d(p)
%
% Input:     p = order of Hermite-chaos
%
% Output:    e = (p+1)x(p+1)x(p+1) matrix containing the result as integers.
% 
% By Dongbin Xiu   5/03/2002
%

np = ceil((3*p+1)/2);

[z,w] = HermiteZW(np);  

factor = sqrt(pi);

J = zeros(p+1, np);
for order=0:p
  J(order+1, :) = HermiteF(z', order);
end
  
for i=1:p+1
  for j=i:p+1
     for k=j:p+1
	  e(i,j,k) = sum(J(i,:).*J(j,:).*J(k,:).*w')/factor;
          e(i,j,k) = round(e(i,j,k));   % e_ijk are integers

	  e(i,k,j) = e(i,j,k);
	  e(j,i,k) = e(i,j,k);
	  e(j,k,i) = e(i,j,k);
	  e(k,i,j) = e(i,j,k);
	  e(k,j,i) = e(i,j,k);
      end
   end
end



