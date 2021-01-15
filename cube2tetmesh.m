function varargout = cube2tetmesh(varargin)
% cube2tetmesh generates uniform tetangle mesh on a cube domain
%
% <SYNTAX>
%   mesh = cube2tetmesh(xl, xr, yl, yr, mx, my)
%   mesh = cube2tetmesh()
%   mesh = cube2tetmesh(mx)
%   mesh = cube2tetmesh(xl, xr)
%   mesh = cube2tetmesh(xl, xr, mx)
%   mesh = cube2tetmesh(xl, xr, yl, yr)
%   mesh = cube2tetmesh(xl, xr, yl, yr, mx)
%   mesh = cube2tetmesh(xl, xr, yl, yr, mx, my)
%
%   [v4tet, x, y] = cube2tetmesh(__)
%
% <DESCRIPTION>
% mesh = cube2tetmesh(xl, xr, yl, yr, mx, my) generates 
% a uniform tetangular mesh on a cubeangle `[xl, xr] x [yl, yr] x [zl x zr]`
% by dividing $x$-, $y$- and $z$-axis by `mx`, `my` and `mz` elements.
% --Example:cube2tetmesh_ex1
%
% mesh = cube2tetmesh() generates `10 x 10 x 10` mesh on the unit cube.
% --Example:cube2tetmesh_ex2
% 
% mesh = cube2tetmesh(mx) generates `mx x mx x mx` mesh on the unit cube.
%
% mesh = cube2tetmesh(xl, xr) generates `10 x 10 x 10` mesh on `[xl, xr]^3`.
%
% mesh = cube2tetmesh(xl, xr, mx) generates `mx x mx x mx` mesh on `[xl, xr]^3`.
%
% mesh = cube2tetmesh(xl, xr, yl, yr, zl, zr) generates `10 x 10 x 10` mesh on `[xl, xr] x [yl, yr] x [zl, zr]`.
%
% mesh = cube2tetmesh(xl, xr, yl, yr, zl, zr, mx) generates `mx x mx x mx` mesh on `[xl, xr] x [yl, yr] x [zl, zr]`.
%
% mesh = cube2tetmesh(xl, xr, yl, yr, zl, zr, mx, my, mz) generates `mx x my x mz` mesh on `[xl, xr] x [yl, yr] x [zl, zr]`.
%
% [v4tet, x, y] = cube2tetmesh(__) returns raw mesh data.
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
%     - zl (double)
%          $z$-min
%     - zr (double)
%          $z$-max
%     - mx (integer)
%          the number of grids on $x$-axis
%     - my (integer)
%          the number of grids on $y$-axis
%     - mz (integer)
%          the number of grids on $z$-axis
%
% <OUTPUT>
%     - mesh (Mesh)
%          Uniform tetrahedral mesh on a cuboid.
%     - v4tet (matetx)
%          Tetrahedron connectivity, specified by 4-row matetx where each column specifies a tetrahedron.
%     - x (vector)
%          $x$-coordinates, specified by a row vector.
%     - y (vector)
%          $y$-coordinates, specified by a row vector.
%     - z (vector)
%          $y$-coordinates, specified by a row vector.
%
% See also cube2cubemesh tetsurfh Mesh

% Copyright 2019 Dohyun Kim / CC BY-NC

% Contact: dhkim.cse@gmail.com
% Developed using MATLAB.ver 9.7 (R2019b) on Microsoft Windows 10 Home

%%
switch nargin
    case 0
        xl = 0; xr = 1; yl = 0; yr = 1; zl = 0; zr = 1;
        mx = 10; my = 10; mz = 10;
    case 1
        mx = varargin{1};
        xl = 0; xr = 1; yl = 0; yr = 1; zl = 0; zr = 1;
        my = mx; mz = mx;
    case 2
        xl = varargin{1}; xr = varargin{2};
        yl = xl; yr = xr;
        zl = xl; zr = xr;
        mx = 10; my = 10;
    case 3
        xl = varargin{1}; xr = varargin{2}; mx = varargin{3};
        yl = xl; yr = xr;
        zl = xl; zr = xr;
        my = mx; mz = mx;
    case 6
        xl = varargin{1}; xr = varargin{2};
        yl = varargin{3}; yr = varargin{4};
        zl = varargin{5}; zr = varargin{6};
        mx = 10; my = 10; mz = 10;
    case 7
        xl = varargin{1}; xr = varargin{2};
        yl = varargin{3}; yr = varargin{4};
        zl = varargin{5}; zr = varargin{6};
        mx = varargin{7}; my = mx; mz = mx;
    case 9
        xl = varargin{1}; xr = varargin{2};
        yl = varargin{3}; yr = varargin{4};
        zl = varargin{5}; zr = varargin{6};
        mx = varargin{7}; my = varargin{8}; mz = varargin{9};
    otherwise
        error('Incorrect number of inputs');
end
            
[x,y,z] = ndgrid(linspace(xl,xr,mx+1), linspace(yl,yr,my+1), linspace(zl,zr,mz+1));
x = x(:); y = y(:); z = z(:);
v4cube = [1, 2, (mx+1) + 1, (mx+1) + 2];
v4cube = v4cube(:) + [0, (mx+1)*(my+1)];
v4cube = v4cube(:);
v4tet = v4cube([1,2,4,6;
                1,4,3,5;
                3,7,5,8;
                4,5,6,8;
                1,4,5,6;
                3,5,4,8].');
v4tet = v4tet(:) + (0:mx-1);
v4tet = v4tet(:) + (0:(mx+1):(mx+1)*(my-1));
v4tet = v4tet(:) + (0:(mx+1)*(my+1):(mx+1)*(my+1)*(mz-1));
v4tet = reshape(v4tet, 4, []);

if nargout < 2
    varargout = {Mesh(v4tet, x, y, z)};
else
    varargout = {v4tet, x, y, z};
end