function z = HermiteZeros(degree)
%
% HermiteZeros.m - Evaluates the zeros of Hermite polynomial
%                  with degree 'degree'.
%
% Syntax :    z = HermiteZeros( degree )
%
% Input  :    degree - degree of Hermite polynomial (number of zeros)
%
% Output :    z - zeros as array (degree x 1)
%
% By Dongbin Xiu   5/03/2002
%

maxit = 50;
EPS = 1.0e-14;
dth = pi/(degree);   % ??? from C code

rlast=0;

% If the degree of the polynomial is zero (or less), no roots
if degree <= 0 
  return;
else
  z = zeros(degree, 1);
end
  
for k = 1:degree
  r = -cos((2*(k-1) + 1) * dth);
  if k ~= 1 
    r = (r + rlast)/2;
  end
  
  for j = 1:maxit
    poly = HermiteF(r, degree);
    pder = HermiteD(r, degree);
      
    delr = -poly./(pder - sum(1./(r-z(1:k-1))).*poly);
    r    = r + delr;
    if  abs(delr) < EPS 
      break;
    end
  end

  if j == maxit 
    fprintf('Maximum number of iteration reached in function: HermiteZeros\n');
  end

  z(k)  = r;
  rlast = r;
end
