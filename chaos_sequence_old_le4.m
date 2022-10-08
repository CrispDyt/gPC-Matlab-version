function q = chaos_sequence(ndim, p)
%
%

if ndim <= 0
   fprintf('Invalid dimensionality! (need n>0)\n');
end

P = nchoosek(ndim+p, p);
%q = zeros(P, ndim);

switch ndim
   case 1
      q = chaos_sequence_1d(ndim, p);
   case 2
      q = chaos_sequence_2d(ndim, p);
   case 3
      q = chaos_sequence_3d(ndim, p);
   case 4
      q = chaos_sequence_4d(ndim, p);
   otherwise
      fprintf('Invalid number of dimensions in chaos_sequence.\n');
end 

% -------------------------------------------------------
% 1d code
% -------------------------------------------------------
function q = chaos_sequence_1d(ndim, p)

P = nchoosek(ndim+p, p);
q = zeros(P,ndim);

q(:,1) = (0:p)';

% -------------------------------------------------------
% 2d code
% -------------------------------------------------------
function q = chaos_sequence_2d(ndim, p)

P = nchoosek(ndim+p, p);
q = zeros(P,ndim);

d=1;
for order=0:p
 for i=order:-1:0
  q(d,1)=i;
  q(d,2)=order-i;
  d=d+1;
 end
end

if (d-1) ~= P
fprintf('Suspicious behavior in chaos_sequence\n');
end

% -------------------------------------------------------
% 3d code
% -------------------------------------------------------
function q = chaos_sequence_3d(ndim, p)

P = nchoosek(ndim+p, p);
q = zeros(P,ndim);

d=1;
for order=0:p
 for i=order:-1:0
  for j=order-i:-1:0
   q(d,1)=i;
   q(d,2)=j;
   q(d,3)=order-i-j;
   d=d+1;
  end
 end
end

if (d-1) ~= P
fprintf('Suspicious behavior in chaos_sequence\n');
end


% -------------------------------------------------------
% 4d code
% -------------------------------------------------------
function q = chaos_sequence_4d(ndim, p)

P = nchoosek(ndim+p, p);
q = zeros(P,ndim);

d=1;
for order=0:p
 for i=order:-1:0
  for j=order-i:-1:0
   for k=order-i-j:-1:0
    q(d,1)=i;
    q(d,2)=j;
    q(d,3)=k;
    q(d,4)=order-i-j-k;
    d=d+1;
   end
  end
 end
end

if (d-1) ~= P
fprintf('Suspicious behavior in chaos_sequence\n');
end
