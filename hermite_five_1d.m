function e = hermite_five_1d(p)
%
% hermite_five_1d.m - Evaluate the inner product of 1d Hermite-chaos with
%                     entries
%
% Syntax     e = hermite_five_1d(p)
%
% Input:     p = order of Hermite-chaos
%
% Output:    e = (p+1)^5 matrix containing the result 
%                as integers.
% 
% By Dongbin Xiu   9/11/2002
%

e = zeros(p+1,p+1,p+1,p+1, p+1);

np = ceil((5*p+1)/2);

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
           for m=l:p+1
	  e(i,j,k,l,m)=sum(J(i,:).*J(j,:).*J(k,:).*J(l,:).*J(m,:).*w')/factor;
          e(i,j,k,l,m) = round(e(i,j,k,l,m));   % e_ijk are integers

	   e(i,k,j,l,m) = e(i,j,k,l,m);
	   e(j,i,k,l,m) = e(i,j,k,l,m);
	   e(j,k,i,l,m) = e(i,j,k,l,m);
	   e(k,i,j,l,m) = e(i,j,k,l,m);
	   e(k,j,i,l,m) = e(i,j,k,l,m);

	   e(i,k,l,j,m) = e(i,j,k,l,m);
	   e(j,i,l,k,m) = e(i,j,k,l,m);
	   e(j,k,l,i,m) = e(i,j,k,l,m);
	   e(k,i,l,j,m) = e(i,j,k,l,m);
	   e(k,j,l,i,m) = e(i,j,k,l,m);

	   e(i,l,k,j,m) = e(i,j,k,l,m);
	   e(j,l,i,k,m) = e(i,j,k,l,m);
	   e(j,l,k,i,m) = e(i,j,k,l,m);
	   e(k,l,i,j,m) = e(i,j,k,l,m);
	   e(k,l,j,i,m) = e(i,j,k,l,m);

	   e(l,i,k,j,m) = e(i,j,k,l,m);
	   e(l,j,i,k,m) = e(i,j,k,l,m);
	   e(l,j,k,i,m) = e(i,j,k,l,m);
	   e(l,k,i,j,m) = e(i,j,k,l,m);
	   e(l,k,j,i,m) = e(i,j,k,l,m);

% ------------------------------------------
	   e(i,k,j,m,l) = e(i,j,k,l,m);
	   e(j,i,k,m,l) = e(i,j,k,l,m);
	   e(j,k,i,m,l) = e(i,j,k,l,m);
	   e(k,i,j,m,l) = e(i,j,k,l,m);
	   e(k,j,i,m,l) = e(i,j,k,l,m);

	   e(i,k,l,m,j) = e(i,j,k,l,m);
	   e(j,i,l,m,k) = e(i,j,k,l,m);
	   e(j,k,l,m,i) = e(i,j,k,l,m);
	   e(k,i,l,m,j) = e(i,j,k,l,m);
	   e(k,j,l,m,i) = e(i,j,k,l,m);

	   e(i,l,k,m,j) = e(i,j,k,l,m);
	   e(j,l,i,m,k) = e(i,j,k,l,m);
	   e(j,l,k,m,i) = e(i,j,k,l,m);
	   e(k,l,i,m,j) = e(i,j,k,l,m);
	   e(k,l,j,m,i) = e(i,j,k,l,m);

	   e(l,i,k,m,j) = e(i,j,k,l,m);
	   e(l,j,i,m,k) = e(i,j,k,l,m);
	   e(l,j,k,m,i) = e(i,j,k,l,m);
	   e(l,k,i,m,j) = e(i,j,k,l,m);
	   e(l,k,j,m,i) = e(i,j,k,l,m);

% ------------------------------------------
	   e(i,k,m,j,l) = e(i,j,k,l,m);
	   e(j,i,m,k,l) = e(i,j,k,l,m);
	   e(j,k,m,i,l) = e(i,j,k,l,m);
	   e(k,i,m,j,l) = e(i,j,k,l,m);
	   e(k,j,m,i,l) = e(i,j,k,l,m);

	   e(i,k,m,l,j) = e(i,j,k,l,m);
	   e(j,i,m,l,k) = e(i,j,k,l,m);
	   e(j,k,m,l,i) = e(i,j,k,l,m);
	   e(k,i,m,l,j) = e(i,j,k,l,m);
	   e(k,j,m,l,i) = e(i,j,k,l,m);

	   e(i,l,m,k,j) = e(i,j,k,l,m);
	   e(j,l,m,i,k) = e(i,j,k,l,m);
	   e(j,l,m,k,i) = e(i,j,k,l,m);
	   e(k,l,m,i,j) = e(i,j,k,l,m);
	   e(k,l,m,j,i) = e(i,j,k,l,m);

	   e(l,i,m,k,j) = e(i,j,k,l,m);
	   e(l,j,m,i,k) = e(i,j,k,l,m);
	   e(l,j,m,k,i) = e(i,j,k,l,m);
	   e(l,k,m,i,j) = e(i,j,k,l,m);
	   e(l,k,m,j,i) = e(i,j,k,l,m);

% ------------------------------------------
	   e(i,m,k,j,l) = e(i,j,k,l,m);
	   e(j,m,i,k,l) = e(i,j,k,l,m);
	   e(j,m,k,i,l) = e(i,j,k,l,m);
	   e(k,m,i,j,l) = e(i,j,k,l,m);
	   e(k,m,j,i,l) = e(i,j,k,l,m);

	   e(i,m,k,l,j) = e(i,j,k,l,m);
	   e(j,m,i,l,k) = e(i,j,k,l,m);
	   e(j,m,k,l,i) = e(i,j,k,l,m);
	   e(k,m,i,l,j) = e(i,j,k,l,m);
	   e(k,m,j,l,i) = e(i,j,k,l,m);

	   e(i,m,l,k,j) = e(i,j,k,l,m);
	   e(j,m,l,i,k) = e(i,j,k,l,m);
	   e(j,m,l,k,i) = e(i,j,k,l,m);
	   e(k,m,l,i,j) = e(i,j,k,l,m);
	   e(k,m,l,j,i) = e(i,j,k,l,m);

	   e(l,m,i,k,j) = e(i,j,k,l,m);
	   e(l,m,j,i,k) = e(i,j,k,l,m);
	   e(l,m,j,k,i) = e(i,j,k,l,m);
	   e(l,m,k,i,j) = e(i,j,k,l,m);
	   e(l,m,k,j,i) = e(i,j,k,l,m);

% ------------------------------------------
	   e(m,i,k,j,l) = e(i,j,k,l,m);
	   e(m,j,i,k,l) = e(i,j,k,l,m);
	   e(m,j,k,i,l) = e(i,j,k,l,m);
	   e(m,k,i,j,l) = e(i,j,k,l,m);
	   e(m,k,j,i,l) = e(i,j,k,l,m);

	   e(m,i,k,l,j) = e(i,j,k,l,m);
	   e(m,j,i,l,k) = e(i,j,k,l,m);
	   e(m,j,k,l,i) = e(i,j,k,l,m);
	   e(m,k,i,l,j) = e(i,j,k,l,m);
	   e(m,k,j,l,i) = e(i,j,k,l,m);

	   e(m,i,l,k,j) = e(i,j,k,l,m);
	   e(m,j,l,i,k) = e(i,j,k,l,m);
	   e(m,j,l,k,i) = e(i,j,k,l,m);
	   e(m,k,l,i,j) = e(i,j,k,l,m);
	   e(m,k,l,j,i) = e(i,j,k,l,m);

	   e(m,l,i,k,j) = e(i,j,k,l,m);
	   e(m,l,j,i,k) = e(i,j,k,l,m);
	   e(m,l,j,k,i) = e(i,j,k,l,m);
	   e(m,l,k,i,j) = e(i,j,k,l,m);
	   e(m,l,k,j,i) = e(i,j,k,l,m);
            end
         end
      end
   end
end



