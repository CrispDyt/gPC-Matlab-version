% CHEBFFT  Chebyshev differentiation via FFT. Simple, not optimal.
%          If v is complex, delete "real" commands.

  function w = chebfft(v)
  N = length(v)-1; if N==0 w=0; break, end
  x = cos((0:N)'*pi/N);
  v = v(:); V = [v; flipud(v(2:N))];      % transform x -> theta
  U = real(fft(V));
  W = real(ifft(1i*[0:N-1 0 1-N:-1]'.*U));
  w = zeros(N+1,1);
  w(2:N) = -W(2:N)./sqrt(1-x(2:N).^2);    % transform theta -> x
  w(1) = sum((0:N-1)'.^2.*U(1:N))/N + .5*N*U(N+1);
  w(N+1) = sum((-1).^(1:N)'.*(0:N-1)'.^2.*U(1:N))/N + ...
             .5*(-1)^(N+1)*N*U(N+1);

