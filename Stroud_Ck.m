function [z,w] = Stroud_Ck( n, m )
%
% This routine returns the Stroud quadrature points for Cn, n-cube.
% Such quadratures are exact for n-dimensional multiple integral 
% in [-1,1]^n with polynomial integrants of degree up to P^m.
% 
% Syntax:  [z,w] = Stroud_Ck( n, m )
%
% Input :  n, dimensionality of the space;
%          m, polynomial exactness.
%
% Only m = 2 and 3 are defined, with number of points n+1 (m=2)
% and 2*n (m=3).  
% 
% Return: z = [npt, n] are the n-dim coordinates for npt points;
%         w = [npt, 1] are the weights.
%         For m=2,3, weights are equal, i.e. w=Volume/npt=2^n/npt.
%
% Formulas are from A.H. Stroud, "Remards on the Disposition of
% Points in Numerical Integration Formulas", Mathematics Tables and
% Other Adis to Computations, Vol. 11(60), 257-261, 1957.
%
% Written and tested by Dongbin Xiu, 12/22/2003.
%

switch m    
 case 2
  npt=n+1;
  z = zeros(npt,n); w=2^n/npt*ones(npt,1);
  
  rr = floor(n/2);
  for k=0:n
    for r=1:rr
      z(k+1, 2*r-1) = sqrt(2/3)*cos(2*r*k*pi/(n+1));
      z(k+1, 2*r)   = sqrt(2/3)*sin(2*r*k*pi/(n+1));
    end
    if mod(n,2) == 1
      z(k+1, n) = (-1)^k/sqrt(3);
    end
  end
 case 3
  npt=2*n;
  z = zeros(npt,n); w=2^n/npt*ones(npt,1);

   rr = floor(n/2);
   for k=1:npt
     for r=1:rr
       z(k, 2*r-1) = sqrt(2/3)*cos((2*r-1)*k*pi/n);
       z(k, 2*r)   = sqrt(2/3)*sin((2*r-1)*k*pi/n);
     end
     if mod(n,2) == 1
       z(k, n) = (-1)^k/sqrt(3);
     end
   end
 otherwise
  fprintf('Stroud_Ck: unknown order of k.\n');
end


