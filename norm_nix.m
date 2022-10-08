function m = norm_nix( k )

%d=length(k);
%m=1;

%for i=1:d
%    %m = m * max(k(i),1);
%    m = m * k(i) + 1;
%end

m = prod( k+1 ) - 1;