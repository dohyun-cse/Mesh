function varargout = line2linemesh(varargin)
% line2linemesh generates uniform line mesh on a interval
%
% <SYNTAX>
%   mesh = line2linemesh(xl, xr, mx)
%   mesh = line2linemesh()
%   mesh = line2linemesh(mx)
%   mesh = line2linemesh(xl, xr)
%
%   [v4e, x] = line2linemesh(__)
%
% <DESCRIPTION>
% mesh = line2linemesh(xl, xr, mx) generates a uniform line mesh 
% on an interval `[xl, xr]` by dividing $x$-axis by `mx`.
%
% mesh = line2linemesh() generates `10` mesh on `[0, 1]`.
% --Example:line2linemesh_ex1
% 
% mesh = line2linemesh(mx) generates `mx` mesh on `[0, 1]`.
%
% mesh = line2linemesh(xl, xr) generates `10` mesh on `[xl, xr]`.
%
% [v4e, x] = line2linemesh(__) returns raw mesh data.
%
% <INPUT>
%     - xl (double)
%          $x$-min
%     - xr (double)
%          $x$-max
%     - mx (integer)
%          the number of grids on $x$-axis
%
% <OUTPUT>
%     - mesh (Mesh)
%          Uniform triangular mesh on an interval.
%     - v4e (matrix)
%          Line connectivity, specified by 2-row matrix where each column specifies a line.
%     - x (vector)
%          $x$-coordinates, specified by a row vector.
%
% See also rect2rectmesh trisurfh Mesh

% Copyright 2019 Dohyun Kim / CC BY-NC

% Contact: dhkim.cse@gmail.com
% Developed using MATLAB.ver 9.7 (R2019b) on Microsoft Windows 10 Home

%%
switch nargin
    case 0
        mx = 10;
        xl = 0;
        xr = 1;
    case 1
        mx = varargin{1};
        xl = 0;
        xr = 1;
    case 2
        xl = varargin{1}; xr = varargin{2};
        mx = 10;
    case 3
        xl = varargin{1}; xr = varargin{2};
        mx = varargin{3};
end

x = linspace(xl, xr, mx+1);

% Left bottom two triangles
v4e = [1:length(x)-1;2:length(x)];
if nargout < 2
    varargout = {Mesh(v4e, x)};
else
    varargout = {v4e, x};
end