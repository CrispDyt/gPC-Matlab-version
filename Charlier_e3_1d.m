function e3 = Charlier_e3_1d(n, alpha)
%
% Charlier_ijk_1d.m - Evaluate the triple inner product of 1d 
%                     Charlier-chaos
%
% Syntax     e3 = Charlier_e3_1d(n, alpha)
%
% Input:     n = order of Charlier-chaos, with parameter alpha.
%
% Output:    e3 = (n+1)x(n+1)x(n+1) matrix containing the result.
% 
% By Dongbin Xiu   5/03/2003
%

np = ceil((3*n+1)/2);

[z,w] = CharlierZW(np, alpha);  

J = zeros(n+1, np);
for i=0:n
  J(i+1, :) = CharlierF(z', i, alpha);
end
  
for i=1:n+1
  for j=i:n+1
     for k=j:n+1
	  e3(i,j,k) = sum(J(i,:).*J(j,:).*J(k,:).*w');
          if(abs(e3(i,j,k)) < 10^(-10))  
             e3(i,j,k)=0;
          end

	  e3(i,k,j) = e3(i,j,k);
	  e3(j,i,k) = e3(i,j,k);
	  e3(j,k,i) = e3(i,j,k);
	  e3(k,i,j) = e3(i,j,k);
	  e3(k,j,i) = e3(i,j,k);
      end
   end
end



