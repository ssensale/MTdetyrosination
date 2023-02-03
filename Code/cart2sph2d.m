function varargout = cart2sph2d(x, y, z)
if nargin == 1
    y = x(:, 2);
    z = x(:, 3);
    x = x(:, 1);
end

% cartesian to spherical conversion
hxy     = hypot(x, y);
rho     = hypot(hxy, z);
theta   = 90 - atan2(z, hxy) * 180 / pi;
phi     = atan2(y, x) * 180 / pi;

% % convert to degrees and theta to colatitude
% theta   = 90 - rad2deg(theta);
% phi     = rad2deg(phi);

% format output
if nargout <= 1
    varargout{1} = [theta phi rho];
    
elseif nargout == 2
    varargout{1} = theta;
    varargout{2} = phi;
    
else
    varargout{1} = theta;
    varargout{2} = phi;
    varargout{3} = rho;
end
    