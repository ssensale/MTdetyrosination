function trans = localToGlobal3d(varargin)
if nargin == 1
    % all components are bundled in  the first argument
    var     = varargin{1};
    center  = var(1:3);
    theta   = var(4);
    phi     = var(5);
    psi     = 0;
    if length(var) > 5
        psi = var(6);
    end
    
elseif nargin == 4
    % arguments = center, then the 3 angles
    center  = varargin{1};
    theta   = varargin{2};
    phi     = varargin{3};
    psi     = varargin{4};    
    
elseif nargin > 4
    % center is given in 3 arguments, then 3 angles
    center  = [varargin{1} varargin{2} varargin{3}];
    theta   = varargin{4};
    phi     = varargin{5};
    psi     = 0;
    if nargin > 5
        psi = varargin{6};    
    end
end
    
% conversion from degrees to radians
k = pi / 180;

% rotation around normal vector axis
rot1    = createRotationOz(psi * k);

% colatitude
rot2    = createRotationOy(theta * k);

% longitude
rot3    = createRotationOz(phi * k);

% shift center
tr      = createTranslation3d(center);

% create final transform by concatenating transforms
trans   = tr * rot3 * rot2 * rot1;
