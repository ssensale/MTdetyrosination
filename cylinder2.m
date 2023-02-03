function [XX,YY,ZZ] = cylinder2(varargin)
narginchk(0,4);
[cax,args,nargs] = axescheck(varargin{:});
% Extract input arguments
n = 20;
r = [1 1]';
d = [0,0,1];
if nargs > 0, r = args{1}; end
if nargs > 1, d = args{2}; end
if nargs > 2, n = args{3}; end
r = r(:); % Make sure r and d are vectors.
d = d(:);
m = length(r); 
if m==1, r = [r;r]; m = 2; end
theta = (0:n)/n*2*pi;
sintheta = sin(theta); sintheta(n+1) = 0;
% Generate basic cylinder
x = r * cos(theta);
y = r * sintheta;
z = (0:m-1)'/(m-1) * ones(1,n+1);
% Define rotation axis and angle
d0      = [0,0,1];
rotaxis = cross(d0,d);
if norm(rotaxis)==0
    rotaxis = [1,0,0];
end
rotaxis = rotaxis/norm(rotaxis);
angle   = -atan2(norm(cross(d0,d)),dot(d0,d));
% Rotation quaternion
q(1:3,1) = rotaxis*sin(angle/2);
q(4,1)   = cos(angle/2);
% Rotation matrix
Q   = [0, -q(3), q(2);q(3), 0, -q(1);-q(2), q(1), 0];
C   = eye(3)*(q(4)^2-q(1:3)'*q(1:3))+2*q(1:3)*q(1:3)'-2*q(4)*Q;
% Generate cylinder
x1      = reshape(x,1,size(x,1)*size(x,2));
y1      = reshape(y,1,size(y,1)*size(y,2));
z1      = reshape(z,1,size(z,1)*size(z,2));
M       = (C*[x1;y1;z1])';
X1      = M(:,1);
Y1      = M(:,2);
Z1      = M(:,3);
X       = reshape(X1,size(x,1),size(x,2));
Y       = reshape(Y1,size(y,1),size(y,2));
Z       = reshape(Z1,size(z,1),size(z,2));
% Output
if nargout == 0
    cax = newplot(cax);
    surf(X,Y,Z,'parent',cax)
else
    XX = X; YY = Y; ZZ = Z;
end
