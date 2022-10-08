function [z,w] = Xiu_Beta( n, m, alpha, beta )
%
% This routine returns the Xiu quadrature points for (-1,1)^n.
% Such quadratures are exact for n-dimensional multiple integral 
% with polynomial exactness of m=2,3, with beta weights 
%     (1-x)^alapha (1+x)^beta
% 
% Syntax:  [z,w] = Xiu_Beta( n, m, alpha, beta )
%
% Input :  n, dimensionality of the space;
%          m, polynomial exactness;
%          alpha, beta >=0.
%
% Only m = 2 and 3 are defined, with number of points n+1 (m=2)
% and 2*n (m=3).  And m=3 is possible only for alpha=beta
% 
% Return: z = [npt, n] are the n-dim coordinates for npt points;
%         w = [npt, 1] are the weights.
%         For m=2,3, weights are equal, i.e. w=Volume/npt=2^n/npt.
%
% Formulas are from D. Xiu, "Numerical Integration Formulas of Degree
% Two", Applied Numerical Mathematics, 2007.
%
% Written and tested by Dongbin Xiu, 7/24/2007.
%

switch m    
 case 2
    [z,w] = Xiu_Hermite(n, m);
    z = (2*sqrt((alpha+1)*(beta+1)/(alpha+beta+3)) * z - (alpha-beta)) / (alpha+beta+2);
 case 3
     if (alpha == beta)
         [z,w] = Xiu_Hermite(n,m);
         z = z/sqrt(2*alpha+3);
     else
         fprintf('Xiu_Beta: 3rd order does not exist.\n');
     end
 otherwise
  fprintf('Xiu_Hermite: undefined order (must be 2 or 3).\n');
end


