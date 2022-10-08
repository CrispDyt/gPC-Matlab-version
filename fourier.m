% FOURIER - compute x = Fourier grid, D = differentiation matrix
% EVEN number of points used.

  function [D,D2,x] = fourier(N)
  if N==0 
    D=0; x=0; 
    return, 
  end
  if mod(N,2)==1 
    fprintf('fourier.m: even points please!\n'); 
    D=0; x=0; return, 
  end
  h = 2*pi/N; x = h*(1:N)';
  column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
  D = toeplitz(column,column([1 N:-1:2]));

  column = [-pi^2/(3*h^2)-1/6 ...
            -.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
  D2 = toeplitz(column);        %  <- 2nd order diff.