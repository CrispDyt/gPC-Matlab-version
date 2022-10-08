function [z,w] = Xiu_Gamma( n, alpha )
%
% This routine returns the Xiu quadrature points for (0,\infty)^n.
% Such quadratures are exact for n-dimensional multiple integral 
% with polynomial exactness of 2, with gamma weights 
%     exp(-x) x^alpha
% The formulas has (n+1) equal weight.
% 
% Syntax:  [z,w] = Xiu_Gamma( n, alpha)
%
% Input :  n, dimensionality of the space;
%          m, polynomial exactness;
%          alpha >=0.
%
% 
% Return: z = [npt, n] are the n-dim coordinates for npt points;
%         w = [npt, 1] are the weights.
%         Weights are equal, i.e. w=Volume/npt=2^n/(n+1).
%
% Formulas are from D. Xiu, "Numerical Integration Formulas of Degree
% Two", Applied Numerical Mathematics, 2007.
%
% Written and tested by Dongbin Xiu, 7/24/2007.
%

    [z,w] = Xiu_Hermite(n, 2);
    z = -sqrt(alpha+1) * z + (alpha+1);
    

