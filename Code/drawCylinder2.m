function varargout = drawCylinder(cyl, varargin)
if iscell(cyl)
    res = zeros(length(cyl), 1);
    for i = 1:length(cyl)
        res(i) = drawCylinder(cyl{i}, varargin{:});
    end
    
    if nargout > 0
        varargout{1} = res;
    end    
    return;
end

% default values
N = 128;
closed = true;

% check number of discretization steps
if ~isempty(varargin)
    var = varargin{1};
    if isnumeric(var)
        N = var;
        varargin = varargin(2:end);
    end
end

% check if cylinder must be closed or open
if ~isempty(varargin)
    var = varargin{1};
    if ischar(var)
        if strncmpi(var, 'open', 4)
            closed = false;
            varargin = varargin(2:end);
        elseif strncmpi(var, 'closed', 5)
            closed = true;
            varargin = varargin(2:end);
        end
    end
end


%% Computation of mesh coordinates

% extreme points of cylinder
p1 = cyl(1:3);
p2 = cyl(4:6);

% radius of cylinder
r = cyl(7);

% compute orientation angle of cylinder
[theta, phi, rho] = cart2sph2d(p2 - p1);
dphi = linspace(0, 2*pi, N+1);

% generate a cylinder oriented upwards
x = repmat(cos(dphi) * r, [2 1]);
y = repmat(sin(dphi) * r, [2 1]);
z = repmat([0 ; rho], [1 length(dphi)]);

% transform points 
trans   = localToGlobal3d(p1, theta, phi, 0);
pts     = transformPoint3d([x(:) y(:) z(:)], trans);

% reshape transformed points
x2 = reshape(pts(:,1), size(x));
y2 = reshape(pts(:,2), size(x));
z2 = reshape(pts(:,3), size(x));


%% Display cylinder mesh

% add default drawing options
varargin = [{'FaceColor', 'b', 'edgeColor', 'none'} varargin];

% plot the cylinder as a surface
    hSurf = surf(x2, y2, z2, varargin{:},'FaceAlpha',1);

% eventually plot the ends of the cylinder
if closed
    ind = find(strcmpi(varargin, 'facecolor'), 1, 'last');
    if isempty(ind)
        color = 'k';
    else
        color = varargin{ind+1};
    end

    patch(x2(1,:)', y2(1,:)', z2(1,:)', color, 'edgeColor', 'none');
    patch(x2(2,:)', y2(2,:)', z2(2,:)', color, 'edgeColor', 'none');
end

% format ouptut
if nargout == 1
    varargout{1} = hSurf;
end