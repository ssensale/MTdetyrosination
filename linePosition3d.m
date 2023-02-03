function pos = linePosition3d(point, line)
dp = bsxfun(@minus, point, line(:,1:3));

% direction vector of the line
vl = line(:, 4:6);

% precompute and check validity of denominator
denom = sum(vl.^2, 2);
invalidLine = denom < eps;
denom(invalidLine) = 1;

% compute position using dot product normalized with norm of line vector.
pos = bsxfun(@rdivide, sum(bsxfun(@times, dp, vl), 2), denom);

% position on a degenerated line is set to 0
pos(invalidLine) = 0;
