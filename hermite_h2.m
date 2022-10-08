function e = hermite_h2( im, jm )
%
% hermite_h2.m - Evaluate the inner product of Hermite-chaos basis
%                <H_im H_jm>, where im,jm are multi-dimensional indices
%                containing its order in each dimension.
%                The exact formula is employed.
%
% Syntax     e = hermite_h2( im, jm )
%
% Input:     im,jm = (1 x M) arrays with each entry denoting its order in
%                    that dimension.
%
% Output:    e = the result in integer form.
% 
% By Dongbin Xiu   2/24/2004.
%

M = length(im);

e = 1;
for m = 1:M
  i = im(m);  j = jm(m); 
  if i ~= j
    e = 0;
    break;
  else
    e = e * factorial(i);
  end
end
