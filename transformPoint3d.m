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

