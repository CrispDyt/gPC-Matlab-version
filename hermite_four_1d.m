function e = hermite_four_1d(p)
%
% hermite_four_1d.m - Evaluate the inner product of 1d Hermite-chaos quadruples
%
% Syntax     e = hermite_four_1d(p)
%
% Input:     p = order of Hermite-chaos
%
% Output:    e = (p+1)x(p+1)x(p+1)x(p+1) matrix containing the result 
%                as integers.
% 
% By Dongbin Xiu   9/11/2002
%

e = zeros(p+1,p+1,p+1,p+1);

np = ceil((4*p+1)/2);

[z,w] = HermiteZW(np);  

factor = sqrt(pi);

J = zeros(p+1, np);
for order=0:p
  J(order+1, :) = HermiteF(z', order);
end
  
for i=1:p+1
  for j=i:p+1
     for k=j:p+1
        for l=k:p+1
	   e(i,j,k,l) = sum(J(i,:).*J(j,:).*J(k,:).*J(l,:).*w')/factor;
           e(i,j,k,l) = round(e(i,j,k,l));   % e_ijk are integers

	   e(i,k,j,l) = e(i,j,k,l);
	   e(j,i,k,l) = e(i,j,k,l);
	   e(j,k,i,l) = e(i,j,k,l);
	   e(k,i,j,l) = e(i,j,k,l);
	   e(k,j,i,l) = e(i,j,k,l);

	   e(i,k,l,j) = e(i,j,k,l);
	   e(j,i,l,k) = e(i,j,k,l);
	   e(j,k,l,i) = e(i,j,k,l);
	   e(k,i,l,j) = e(i,j,k,l);
	   e(k,j,l,i) = e(i,j,k,l);

	   e(i,l,k,j) = e(i,j,k,l);
	   e(j,l,i,k) = e(i,j,k,l);
	   e(j,l,k,i) = e(i,j,k,l);
	   e(k,l,i,j) = e(i,j,k,l);
	   e(k,l,j,i) = e(i,j,k,l);

	   e(l,i,k,j) = e(i,j,k,l);
	   e(l,j,i,k) = e(i,j,k,l);
	   e(l,j,k,i) = e(i,j,k,l);
	   e(l,k,i,j) = e(i,j,k,l);
	   e(l,k,j,i) = e(i,j,k,l);
         end
      end
   end
end



