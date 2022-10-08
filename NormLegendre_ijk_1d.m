function e = NormLegendre_ijk_1d(p)
%
% NormLegendre_ijk_1d.m - Evaluate the inner product of 1d normalized
%                   Legendre-triplets.
%
% Syntax     e = NormLegendre_ijk_1d(p)
%
% Input:     p = order of normalized Legendre-chaos
%
% Output:    e = (p+1)x(p+1)x(p+1) matrix containing the result.
% 
% NO WARNING MESSAGE IS GIVEN WHEN PAPAMETERS ARE OUT OF RANGE.
%
% By Dongbin Xiu   06/13/2005
%

tolerance=1e-10;

np = ceil((3*p+1)/2);

alpha=0;  beta = 0;
[z,w] = zwgj(np, alpha, beta);  

factor = 2^(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)/gamma(alpha+beta+2);

J = zeros(p+1, np);
for order=0:p
  J(order+1, :) = LegendreF_Normal(z', order);
end
  
for i=1:p+1
  for j=i:p+1
     for k=j:p+1
	  e(i,j,k) = sum(J(i,:).*J(j,:).*J(k,:).*w')/factor;
          if abs(e(i,j,k)-round(e(i,j,k))) < tolerance
             e(i,j,k) = round(e(i,j,k));
          end

	  e(i,k,j) = e(i,j,k);
	  e(j,i,k) = e(i,j,k);
	  e(j,k,i) = e(i,j,k);
	  e(k,i,j) = e(i,j,k);
	  e(k,j,i) = e(i,j,k);
      end
   end
end



