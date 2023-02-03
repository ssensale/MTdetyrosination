function points = intersectLineCylinder(line, cylinder, varargin)
% default arguments
checkBounds = true;

% parse inputs
while length(varargin)>1
    var = varargin{1};
    if strcmpi(var, 'checkbounds')
        checkBounds = varargin{2};
    else
        error(['Unkown argument: ' var]);
    end
    varargin(1:2) = [];
end


%% Parse cylinder parameters

% Starting point of the line
l0 = line(1:3)';

% Direction vector of the line
dl = line(4:6)';

% Starting position of the cylinder
c0 = cylinder(1:3)';

% Direction vector of the cylinder
dc = cylinder(4:6)' - c0;

% Radius of the cylinder
r = cylinder(7);


%% Resolution of a quadratic equation to find the increment

% Substitution of parameters
e = dl - (dot(dl,dc)/dot(dc,dc))*dc;
f = (l0-c0) - (dot(l0-c0,dc)/dot(dc,dc))*dc;

% Coefficients of 2-nd order equation
A = dot(e, e);
B = 2*dot(e,f);
C = dot(f,f) - r^2;

% compute discriminant
delta = B^2 - 4*A*C;

% check existence of solution(s)
if delta<0
    points = zeros(0, 3);
    return;
end

% extract roots
x1 = (-B + (delta)^0.5)/(2*A);
x2 = (-B - (delta)^0.5)/(2*A);
x = [x1;x2];


%% Estimation of points position

% process the smallest position
x1 = min((x));

% Point on the line: l0 + x*dl = p
point1 = l0 + x1*dl;

% process the greatest position
x2 = max((x));

% Point on the line: l0 + x*dl = p
point2 = l0 + x2*dl;

% Format result
points = [point1' ; point2'];


%% Check if points are located between bounds

if checkBounds
    % cylinder axis
    axis = [c0' dc'];
    
    % compute position on axis
    ts = linePosition3d(points, axis);
    
    % check bounds
    ind = ts>=0 & ts<=1;
    points = points(ind, :);
end