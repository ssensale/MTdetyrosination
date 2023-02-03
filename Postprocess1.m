clc; 
clear;
clf;
tic;
source_dir = pwd; 
d = dir([source_dir, '/Chun*.dat']);
numscreens = length(d);
start_ind = 2;
end_ind=numscreens;
widthgaus=2.35;
diftime=500;
doplots=1;
dohistograms=1;
dofitting=1;
for toma2 = end_ind:-1:start_ind%numscreens
    delta = 0.6;
    Cilindros = importdata('Cylinder.dat');
    Direct = dir('Len*.dat');
    if size(Direct,1) > 1
        Largo = importdata('Lengths.dat');    
    else
        Largo = importdata(Direct.name);    
    end
    Largos = Largo(toma2,:);
    for i = 1:size(Largos,2)
        Largos(i) = min(Largos(i),12);
    end
    multip = 5;
    Detyros = importdata(strcat('Chunks',num2str(diftime * (toma2 - 1) * multip), '.dat'));
    for i = 1:size(Largos,2)
        Detyros(i,2:end) = sort(Detyros(i,2:end));
    end
    radiusin = 0.5;
    radiusout = 12.5;
    discretization = (radiusout - radiusin) * 1000 / delta;
    
    %Place the sites 
    finegrid = delta / 1000;
    pepe = 0;
    for i = 1:size(Largos,2)
       dx = (Cilindros(i,5:7) - Cilindros(i,2:4)) * finegrid / (norm((Cilindros(i,5:7) - Cilindros(i,2:4))));        
       if Largos(i) == 0
           sites{i} = [];
           percdet(i,1:3) = NaN;
           modificados{i} = [];
           newdiameter(i) = 0;
           primero2(i,:) = zeros(1,3);
           ultimo2(i,:) = zeros(1,3);
           largodetyzx(i)=0;
       else
           pepe = pepe+1;
           numerosites(i) = floor(Largos(i)*1000/delta);
           numerodets(i) = nnz(~isnan(Detyros(i,:))) - 1;
           Larger(pepe) = Largos(i);
           numee(pepe) = numerodets(i) * 0.6 / 1000;
           largodetyzx(i)=numerodets(i)*0.6/1000;
           %extend the MT start
           a = floor(lognrnd(2.7952,0.64166,[numerosites(i),1]));  
           comienzo(1) = 1;    
           valores(1:a(1)) = 1;
           for qu = 2:numerosites(i)
                comienzo(qu) = comienzo(qu-1)+a(qu-1);
                valores(comienzo(qu):comienzo(qu)+a(qu)-1) = qu;
           end            
           r = normrnd(0,0.03/widthgaus,[sum(a,1),1]); 
           theta = rand(sum(a,1),1).*2.*pi;
           phi = rand(sum(a,1),1).*pi;
           Px = r.*cos(theta).*sin(phi);
           Py = r.*sin(theta).*sin(phi);
           Pz = r.*cos(phi);
           for qu = 1:sum(a,1)
                proyecto = (Px(qu) * dx(1) + Py(qu) * dx(2) + Pz(qu) * dx(3)) / norm(dx);
                projeccion(qu) = floor(proyecto / (delta / 1000));
                projeccion2(qu) = valores(qu) + floor(proyecto / (delta / 1000));
                resto(qu) = sqrt(r(qu)^2 - proyecto^2);
           end
           primero(i) = min(projeccion2);
           ultimo(i) = max(projeccion2);
           sitioss(i) = ultimo(i) - primero(i) + 1;
           primero2(i,:) = Cilindros(i,2:4) + dx * min(0,primero(i));
           ultimo2(i,:) = primero2(i,:) + dx * max(0,ultimo(i));    
           newdiameter(i) = mean(resto(qu)) + 25 / 1000;            
           for qu = 1:numerodets(i)
                comi = comienzo(Detyros(i,qu+1));
                fini = comi+a(Detyros(i,qu+1))-1;
                proyector{i,qu} = projeccion2(comi:fini);
                proyector{i,qu}= min(proyector{i,qu}):max(proyector{i,qu});%this line is questionable-we either sample the points or we consider interval
           end
           if numerodets(i) > 0
                modificados{i} = -primero(i) + 1 + unique(sort([proyector{i,:}],2),'first');%new detyrosinated sites
           else
                modificados{i} = [];
           end
           clearvars dx sites comienzo valores a r theta phi Px Py Pz qu proyecto projeccion projeccion2 resto
       end
    end
    for i = 1:size(Largos,2)
        modificados{i} = [i,modificados{i}];
        modiphony{i} = [i,primero2(i,1),primero2(i,2),primero2(i,3),ultimo2(i,1),ultimo2(i,2),ultimo2(i,3),norm(ultimo2(i,:)-primero2(i,:)),newdiameter(i)];
    end
   
    clearvars -except fig percdet modificados modiphony start_ind toma2 discretization radiusout delta multip radiusin Largos RATIOS proyector numerodets numerosites sitioss Larger numee doplots dohistograms largodetyzx diftime widthgaus dofitting
    
    Detyros = modificados';

    %go row by row, if any two cells of proyector (i,:) have elements in
    %common, then merge them

    for i=1:size(proyector,1)
        A=proyector(i,:);
        ii = 1;
        while ii<=numel(A)
            kk = [];
            for jj = ii+1:numel(A)
              if and(not(isempty(A{jj})),not(isempty(A{ii})))
                   if not(or(and(A{ii}(1)<A{jj}(1),A{ii}(end)<A{jj}(1)),and(A{ii}(1)>A{jj}(end),A{ii}(end)>A{jj}(end))))
                        BE=unique(cat(2,A{ii},A{jj}));
                        A{ii} = BE;
                        kk = [kk,jj];   
                    end
                end
            end
            A(kk) = []; % remove cells
            ii = ii+1;
        end
        A=A(~cellfun('isempty',A));
        for k=1:size(A,2)
            proyector{i,k}=A{k};
        end
        for k=size(A,2)+2:size(proyector,2)
            proyector{i,k}=[];
        end
    end    

    for i = 1:size(modiphony,2)
        dx(i,:) = (modiphony{i}([5:7]) - modiphony{i}([2:4])) / modiphony{i}(8);
        largositio(i)=0;
        if not(isnan(dx(i,1)))
            if i<=size(proyector,1)
                for j = 1:size(proyector,2)
                    if not(isempty(proyector{i,j}))
                        positionx{i,j}(:,1) = modiphony{i}([2:4]) + dx(i,1:3) .* ...
                            (proyector{i,j}(1)-1) * 0.6 / 1000; %beginning
                        positionx{i,j}(:,2) = modiphony{i}([2:4]) + dx(i,1:3) .* ...
                            proyector{i,j}(end) * 0.6 / 1000; %end
                        largositio(i)=largositio(i)+(proyector{i,j}(end)-proyector{i,j}(1)-1)*0.6/1000;
                    end
                end
            end
        end
    end

    %Print the histograms to file
    for i = 1:size(modiphony,2)
        %Output(i,1)=sitioss(i);%sites available after jittering
        %Output(i,2)=numerosites(i);%sites available before jittering

        Output(i,1)=modiphony{1,i}(8);%MT length after jittering
        Output(i,4)=Largos(i);%MT length before jittering
        Output(i,2)=largositio(i);%length det after jittering
        Output(i,5)=largodetyzx(i);%length det before jittering
        if Largos(i)>0
            Output(i,3)=100*Output(i,2)/Output(i,1);
        else
            Output(i,3)=NaN;
        end
        if Largos(i)>0
            Output(i,6)=100*Output(i,5)/Output(i,4);
        else
            Output(i,6)=NaN;
        end        
    end
    
    OutputMTnozero(:,:)=Output(Output(:,1)>0,:);
    OutputDETnozero(:,:)=Output(Output(:,5)>0,:);

    writematrix(Output,strcat('zHistogramsfor',num2str(toma2),'.dat'),'Delimiter','tab'); 

    %Calculate the net values
    Output2(1)=diftime * (toma2 - 1) * multip/60/60;%time in minutes
    Output2(2)=mean(Output(:,1));%Total MT length after jitt
    Output2(3)=mean(Output(:,4));%Total MT length before jitt
    Output2(4)=sum(Output(:,2));%Total length det after jitt
    Output2(5)=sum(Output(:,5));%Total length det before jitt
    Output2(6)=100*Output2(4)/sum(Output(:,1));%Perc det after jitt
    Output2(7)=100*Output2(5)/sum(Output(:,4));%Perc det before jitt
 
    writematrix(Output2,strcat('zData',num2str(toma2),'.dat'),'Delimiter','tab'); 

 
    clearvars -except source_dir d numscreens start_ind toma2 end_ind doplots dohistograms diftime widthgaus dofitting
    toc
end

function varargout = drawCylinder(cyl, varargin)
%% Input argument processing

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
end
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
end
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
end
function trans = createRotationOz(varargin)
% default values
dx = 0;
dy = 0;
dz = 0;
theta = 0;

% get input values
if length(varargin)==1
    % only angle
    theta = varargin{1};
elseif length(varargin)==2
    % origin point (as array) and angle
    var = varargin{1};
    dx = var(1);
    dy = var(2);
    dz = var(3);
    theta = varargin{2};
elseif length(varargin)==3
    % origin (x and y) and angle
    dx = varargin{1};
    dy = varargin{2};
    dz = varargin{3};
    theta = varargin{3};
end

% compute coefs
cot = cos(theta);
sit = sin(theta);

% create transformation
trans = [...
    cot -sit 0 0;...
    sit  cot 0 0;...
    0 0 1 0;...
    0 0 0 1];

% add the translation part
t = [1 0 0 dx;0 1 0 dy;0 0 1 dz;0 0 0 1];
trans = t*trans/t;
end
function trans = createRotationOy(varargin)
% default values
dx = 0;
dy = 0;
dz = 0;
theta = 0;

% get input values
if length(varargin)==1
    % only angle
    theta = varargin{1};
elseif length(varargin)==2
    % origin point (as array) and angle
    var = varargin{1};
    dx = var(1);
    dy = var(2);
    dz = var(3);
    theta = varargin{2};
elseif length(varargin)==3
    % origin (x and y) and angle
    dx = varargin{1};
    dy = 0;
    dz = varargin{2};
    theta = varargin{3};
elseif length(varargin)==4
    % origin (x and y) and angle
    dx = varargin{1};
    dy = varargin{2};
    dz = varargin{3};
    theta = varargin{4};
end

% compute coefs
cot = cos(theta);
sit = sin(theta);

% create transformation
trans = [...
    cot  0  sit  0;...
    0    1    0  0;...
    -sit 0  cot  0;...
    0    0    0  1];

% add the translation part
t = [1 0 0 dx;0 1 0 dy;0 0 1 dz;0 0 0 1];
trans = t*trans/t;
end
function trans = createTranslation3d(varargin)
if isempty(varargin)
    % assert translation with null vector
    dx = 0;
    dy = 0;
    dz = 0;
elseif length(varargin)==1
    % translation vector given in a single argument
    var = varargin{1};
    dx = var(1);
    dy = var(2);
    dz = var(3);
else
    % translation vector given in 3 arguments
    dx = varargin{1};
    dy = varargin{2};
    dz = varargin{3};
end

% create the translation matrix
trans = [1 0 0 dx ; 0 1 0 dy ; 0 0 1 dz; 0 0 0 1];
end
function varargout = transformPoint3d(varargin)
meshStructSwitch = false;
if length(varargin) == 2
    if isstruct(varargin{1}) && isfield(varargin{1}, 'vertices')
        % If first argument is a struct with the field 'vertices' the 
        % output will be the same struct, but with the transformed vertices
        meshStructSwitch = true;
        meshStruct = varargin{1};
        varargin{1}=varargin{1}.vertices;
    end
    % Point coordinates are given in a single N-by-3-by-M-by-etc argument.
    % Preallocate x, y, and z to size N-by-1-by-M-by-etc, then fill them in
    dim = size(varargin{1});
    dim(2) = 1;
    [x, y, z] = deal(zeros(dim, class(varargin{1})));
    x(:) = varargin{1}(:,1,:);
    y(:) = varargin{1}(:,2,:);
    z(:) = varargin{1}(:,3,:);
    trans  = varargin{2};  
elseif length(varargin) == 4
    % Point coordinates are given in 3 different arrays
    x = varargin{1};
    y = varargin{2};
    z = varargin{3};
    dim = size(x);
    trans = varargin{4};
end

% eventually add null translation
if size(trans, 2) == 3
    trans = [trans zeros(size(trans, 1), 1)];
end

% eventually add normalization
if size(trans, 1) == 3
    trans = [trans;0 0 0 1];
end

% convert coordinates
NP  = numel(x);
try
    % vectorial processing, if there is enough memory
    %res = (trans*[x(:) y(:) z(:) ones(NP, 1)]')';
    %res = [x(:) y(:) z(:) ones(NP, 1)]*trans';    
    res = [x(:) y(:) z(:) ones(NP,1,class(x))] * trans';
    
    % Back-fill x,y,z with new result (saves calling costly reshape())
    x(:) = res(:,1);
    y(:) = res(:,2);
    z(:) = res(:,3);
catch ME
    disp(ME.message)
    % process each point one by one, writing in existing array
    for i = 1:NP
        res = [x(i) y(i) z(i) 1] * trans';
        x(i) = res(1);
        y(i) = res(2);
        z(i) = res(3);
    end
end

% process output arguments
if nargout <= 1
    % results are stored in a unique array
    if length(dim) > 2 && dim(2) > 1
        warning('geom3d:shapeMismatch',...
            'Shape mismatch: Non-vector xyz input should have multiple x,y,z output arguments. Cell {x,y,z} returned instead.')
        varargout{1} = {x,y,z};
    else
        if meshStructSwitch
            meshStruct.vertices = [x y z];
            varargout{1}=meshStruct;
        else
            varargout{1} = [x y z];
        end
    end
    
elseif nargout == 3
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = z;
end
end

