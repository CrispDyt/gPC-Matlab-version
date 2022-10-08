% CHEB - compute x = Chebyshev grid, D = differentiation matrix

  function [D,x] = cheb(N)
  if N==0 D=0; x=1; return, end
  ii = (0:N)';
  x = cos(pi*ii/N);
  c = [2; ones(N-1,1); 2];
  D = zeros(N+1,N+1);
  for j = 0:N;
    denom = c(j+1)*(x(ii+1)-x(j+1)); denom(j+1) = inf;
    D(ii+1,j+1) = c.*(-1).^(ii+j)./denom;
  end
  D = D - diag(sum(D'));

