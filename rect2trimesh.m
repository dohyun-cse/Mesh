function varargout = rect2trimesh(varargin)
% rect2trimesh generates uniform triangle mesh on a square domain
%
% <SYNTAX>
%   mesh = rect2trimesh(xl, xr, yl, yr, mx, my)
%   mesh = rect2trimesh()
%   mesh = rect2trimesh(mx)
%   mesh = rect2trimesh(xl, xr)
%   mesh = rect2trimesh(xl, xr, mx)
%   mesh = rect2trimesh(xl, xr, yl, yr)
%   mesh = rect2trimesh(xl, xr, yl, yr, mx)
%   mesh = rect2trimesh(xl, xr, yl, yr, mx, my)
%
%   [v4e, x, y] = rect2trimesh(__)
%
% <DESCRIPTION>
% mesh = rect2trimesh(xl, xr, yl, yr, mx, my) generates 
% a uniform triangular mesh on a rectangle `[xl, xr] x [yl, yr]` by dividing 
% $x$- and $y$-axis by `mx` and `my` elements.
% --Example:rect2trimesh_ex1
%
% mesh = rect2trimesh() generates `10 x 10` mesh on the unit square.
% --Example:rect2trimesh_ex2
% 
% mesh = rect2trimesh(mx) generates `mx x mx` mesh on the unit square.
%
% mesh = rect2trimesh(xl, xr) generates `10 x 10` mesh on `[xl, xr] x [xl, xr]`.
%
% mesh = rect2trimesh(xl, xr, mx) generates `mx x mx` mesh on `[xl, xr] x [yl, yr]`.
%
% mesh = rect2trimesh(xl, xr, yl, yr) generates `10 x 10` mesh on `[xl, xr] x [yl, yr]`.
%
% mesh = rect2trimesh(xl, xr, yl, yr, mx) generates `mx x mx` mesh on `[xl, xr] x [yl, yr]`.
%
% mesh = rect2trimesh(xl, xr, yl, yr, mx, my) generates `mx x my` mesh on `[xl, xr] x [yl, yr]`.
%
% [v4e, x, y] = rect2trimesh(__) returns raw mesh data.
%
% <INPUT>
%     - xl (double)
%          $x$-min
%     - xr (double)
%          $x$-max
%     - yl (double)
%          $y$-min
%     - yr (double)
%          $y$-max
%     - mx (integer)
%          the number of grids on $x$-axis
%     - my (integer)
%          the number of grids on $y$-axis
%
% <OUTPUT>
%     - mesh (Mesh)
%          Uniform triangular mesh on a rectangle.
%     - v4e (matrix)
%          Triangle connectivity, specified by 3-row matrix where each column specifies a triangle.
%     - x (vector)
%          $x$-coordinates, specified by a row vector.
%     - y (vector)
%          $y$-coordinates, specified by a row vector.
%
% See also rect2rectmesh trisurfh Mesh

% Copyright 2019 Dohyun Kim / CC BY-NC

% Contact: dhkim.cse@gmail.com
% Developed using MATLAB.ver 9.7 (R2019b) on Microsoft Windows 10 Home

%%
switch nargin
    case 0
        xl = 0; xr = 1; yl = 0; yr = 1;
        mx = 10; my = 10;
    case 1
        mx = varargin{1};
        xl = 0; xr = 1; yl = 0; yr = 1;
        my = mx;
    case 2
        xl = varargin{1}; xr = varargin{2};
        yl = xl; yr = xr;
        mx = 10; my = 10;
    case 3
        xl = varargin{1}; xr = varargin{2}; mx = varargin{3};
        yl = xl; yr = xr;
        my = mx;
    case 4
        xl = varargin{1}; xr = varargin{2};
        yl = varargin{3}; yr = varargin{4};
        mx = 10; my = 10;
    case 5
        xl = varargin{1}; xr = varargin{2};
        yl = varargin{3}; yr = varargin{4};
        mx = varargin{5}; my = mx;
    case 6
        xl = varargin{1}; xr = varargin{2};
        yl = varargin{3}; yr = varargin{4};
        mx = varargin{5}; my = varargin{6};
end

[x,y] = ndgrid(linspace(xl, xr, mx+1), linspace(yl, yr, my+1));
x = x(:).'; y = y(:).';

% Left bottom two triangles
v4e = [1, 2, mx+2, mx+3, mx+2, 2];
% Bottom strip
v4e = v4e(:) + (0:mx-1);
% whole domain
v4e = v4e(:) + (0:(mx+1):(mx+1)*(my-1));
% reshape
v4e = reshape(v4e, 3, []);
if nargout < 2
    varargout = {Mesh(v4e, x, y)};
else
    varargout = {v4e, x, y};
end