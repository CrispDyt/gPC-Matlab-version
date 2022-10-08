function poly = CharlierF(x, degree, alpha)
%
% CharlierF.m - Compute the Charlier polynomial C^n(x;alpha) using the
%               backward shift operator.
%
% Syntax:   poly = CharlierF(x, degree, alpha);
%
% Input :   x = x-coordinate in matrix form
%           degree = order of polynomial
%           alpha > 0 is the parameter
%
% Output:   poly = polynomial values at location x stored in same form as x.
%
% by Dongbin Xiu   5/02/2003
%

switch degree
  case -1
    poly = zeros(size(x));
  case 0
    poly = ones(size(x));
  otherwise
    poly = x.*CharlierF(x-1, degree-1, alpha) - ...
           alpha*CharlierF(x, degree-1, alpha);
end
