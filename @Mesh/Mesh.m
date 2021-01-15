classdef Mesh < handle & matlab.mixin.CustomDisplay
% Mesh is a class for n-dimensional simplex mesh
%
% <SYNTAX>
%   mesh = Mesh(v4line, x)
%   mesh = Mesh(v4e, x, y)
%   mesh = Mesh(v4tet, x, y, z)
%
% <DESCRIPTION>
%   mesh = Mesh(v4line, x) generates 1-Dimensional linear mesh information.
%
%   mesh = Mesh(v4e, x, y) generates 2-Dimensional triangular mesh
%   information.
%
%   mesh = Mesh(v4tet, x, y, z) generates 3-Dimensional tetrahedral mesh
%   information.
%
% <INPUT>
%   - v4e (`(n+1) x nrE` matrix)
%       Vertex index for elements, specified by a matrix, where `v4e(:,i)` is
%       the vertex indices of `i`-th element.
%   - x (`1 x nrV` vector)
%       $x$-coordinates of each vertex, specified by a vector.
%   - y (`1 x nrV` vector)
%       $y$-coordinates of each vertex, specified by a vector.
%   - z (`1 x nrV` vector)
%       $z$-coordinates of each vertex, specified by a vector.
%          
% <OUTPUT>
%   - mesh (Mesh)
%       Resulting mesh class
%
% See also line2linemesh rect2trimesh cube2tetmesh

% Copyright 2019 Dohyun Kim / CC BY-NC

% Contact: dhkim.cse@gmail.com
% Developed using MATLAB.ver 9.7 (R2019b) on Microsoft Windows 10 Enterprise

%%
    properties
        
        dim
        shape
        v4e
        x, y, z
        
        
        J
        rx, sx, mx, ry, sy, my, rz, sz, mz
        xr, yr, zr, xs, ys, zs, xm, ym, zm
        
        
        fmask, fmask_neg
        
        v4fb, e4fb, ef4fb
        v4f0, e4f0, ef4f0
        
        J4f0, J4fb
        nx4fb, ny4fb, nz4fb
        nx4f0, ny4f0, nz4f0
        
        BC

        gpuFlag = false

        indexMesh
        
        isFacetBuilt = false
        isTriAnalyzed = false
        isFacetAnalyzed = false
        
    end
        
    properties (Dependent)
        
        x4e, y4e, z4e
        
        x4fb, y4fb, z4fb
        x4f0, y4f0, z4f0
        x4f, y4f, z4f
        
        v4f
        e4f, ef4f
        
        J4f
        nx4f, ny4f, nz4f
        tx4fb, ty4fb
        tx4f0, ty4f0
        
        f4e, e4e
        
        nrE, nrV, nrfb, nrf0, nrf
        
    end
    
    
    
    %%
    methods
        
        function mesh = Mesh(v4e, varargin)
            % Constructor of mesh
            if ~nargin
                return;
            end
            global GPU_ON
            mesh.v4e = v4e;
            mesh.dim = length(varargin);
            ax = 'xyz';
            if GPU_ON
                for i = 1 : length(varargin)
                    mesh.(ax(i)) = gpuArray(varargin{i});
                end
                mesh.gpuFlag = true;
            else
                for i = 1 : length(varargin)
                    mesh.(ax(i)) = varargin{i};
                end
            end
            if mesh.dim + 1 ~= size(mesh.v4e, 1)
                error('Mesh class only takes simplex')
            end
            
            switch mesh.dim
                case 1
                    mesh.fmask = [1,2];
                    mesh.fmask_neg = [1,2];
                    mesh.shape = 'line';
                case 2
                    mesh.fmask = [1,2;2,3;3,1].';
                    mesh.fmask_neg = flipud(mesh.fmask);
                    mesh.shape = 'triangle';
                case 3
                    mesh.fmask = [3,2,1;1,2,4;2,3,4;3,1,4].';
                    mesh.fmask_neg = [1,3,2;2,1,3;3,2,1].';
                    mesh.fmask_neg = [mesh.fmask(mesh.fmask_neg(:,1),:), ...
                                      mesh.fmask(mesh.fmask_neg(:,2),:), ...
                                      mesh.fmask(mesh.fmask_neg(:,3),:)];
                    mesh.shape = 'tetrahedron';
            end
            
            mesh.BC = containers.Map('KeyType', 'char', 'ValueType', 'any');
            mesh.indexMesh = mesh;
        end
        
        h = show(obj, varargin)
        
        function setBC(mesh, bcName, bc)
            if ~ischar(bcName)
                error('Boundary condition name should be char');
            end
            if isa(bc, 'logical')
                mesh.BC(bcName) = bc;
            elseif isa(bc, 'function_handle')
                direction = 'xyz';
                args = arrayfun(@(d) mesh.([d '4fb']), direction(1:mesh.dim), 'un', 0);
                mesh.BC(bcName) = bc(args{:});
                if ~islogical(mesh.BC(bcName))
                    error('function handle for boundary condition should return logical value');
                end
            else
                error('Boundary condition should be provided as logical or function_handle');
            end
        end
        
        function varargout = get(mesh, varargin)
            if length(varargin) ~= nargout
                error('The number of requested properties and the number of output arguments should be the same.')
            end
            varargout = cell(size(varargin));
            for i = 1 : length(varargin)
                varargout{i} = mesh.(varargin{i});
            end
        end
        
        function resetCoordinate(mesh)
            mesh.isTriAnalyzed = false;
            mesh.isFacetAnalyzed = false;
        end
        function resetConnectivity(mesh)
            mesh.isTriAnalyzed = false;
            mesh.isFacetAnalyzed = false;
            mesh.isFacetBuilt = false;
        end
        
        %% SAFE GET METHODS
        
    
        function set.x(mesh, x), mesh.x = x(:).'; resetCoordinate(mesh); end
        function set.y(mesh, y), mesh.y = y(:).'; resetCoordinate(mesh); end
        function set.z(mesh, z), mesh.z = z(:).'; resetCoordinate(mesh); end
        function set.v4e(mesh, v4e), mesh.v4e = v4e; resetConnectivity(mesh); end
        
        % Triangle affine information
        function J  = get.J(mesh),  if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, J  = mesh.J;  end
        function rx = get.rx(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, rx = mesh.rx; end
        function ry = get.ry(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, ry = mesh.ry; end
        function rz = get.rz(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, rz = mesh.rz; end
        function sx = get.sx(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, sx = mesh.sx; end
        function sy = get.sy(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, sy = mesh.sy; end
        function sz = get.sz(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, sz = mesh.sz; end
        function mx = get.mx(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, mx = mesh.mx; end
        function my = get.my(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, my = mesh.my; end
        function mz = get.mz(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, mz = mesh.mz; end
        
        function xr = get.xr(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, xr = mesh.xr; end
        function yr = get.yr(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, yr = mesh.yr; end
        function zr = get.zr(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, zr = mesh.zr; end
        function xs = get.xs(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, xs = mesh.xs; end
        function ys = get.ys(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, ys = mesh.ys; end
        function zs = get.zs(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, zs = mesh.zs; end
        function xm = get.xm(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, xm = mesh.xm; end
        function ym = get.ym(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, ym = mesh.ym; end
        function zm = get.zm(mesh), if ~mesh.isTriAnalyzed, mesh.affineSmplx(); end, zm = mesh.zm; end
        
        % The number
        function nrE = get.nrE(mesh), nrE = size(mesh.v4e, 2); end
        function nrV = get.nrV(mesh), nrV = length(mesh.x); end
        function nrf  = get.nrf(mesh), nrf  = size(mesh.v4fb, 2) + size(mesh.v4f0, 2);  end
        function nrfb = get.nrfb(mesh), nrfb = size(mesh.v4fb, 2); end
        function nrf0 = get.nrf0(mesh), nrf0 = size(mesh.v4f0, 2); end
        
        % facet
        function v4fb = get.v4fb(mesh), if ~mesh.isFacetBuilt, mesh.buildFacets(); end, v4fb = mesh.v4fb; end
        function v4f0 = get.v4f0(mesh), if ~mesh.isFacetBuilt, mesh.buildFacets(); end, v4f0 = mesh.v4f0; end
        function v4f = get.v4f(mesh), v4f = [mesh.v4fb, mesh.v4f0]; end
        
        % Connectivity information
        function e4fb  = get.e4fb(mesh),  if ~mesh.isFacetBuilt, mesh.buildFacets(); end, e4fb  = mesh.e4fb;  end
        function ef4fb = get.ef4fb(mesh), if ~mesh.isFacetBuilt, mesh.buildFacets(); end, ef4fb = mesh.ef4fb; end
        function e4f0  = get.e4f0(mesh),  if ~mesh.isFacetBuilt, mesh.buildFacets(); end, e4f0  = mesh.e4f0;  end
        function ef4f0 = get.ef4f0(mesh), if ~mesh.isFacetBuilt, mesh.buildFacets(); end, ef4f0 = mesh.ef4f0; end
        
        function e4f = get.e4f(mesh), e4f = [repmat(mesh.e4fb,2,1), mesh.e4f0]; end
        function ef4f = get.ef4f(mesh), ef4f = [repmat(mesh.ef4fb,2,1), mesh.ef4f0]; end
        
        function f4e = get.f4e(mesh)
            if ~mesh.isFacetBuilt, mesh.buildFacets(); end
            nrF4e = size(mesh.fmask,2);
            f4e = zeros(nrF4e, mesh.nrE);
            f4e((mesh.e4fb-1)*nrF4e + mesh.ef4fb) = 1:mesh.nrfb;
            f4e((mesh.e4f0-1)*nrF4e + mod(mesh.ef4f0-1, nrF4e)+1) = repmat(mesh.nrfb+1:mesh.nrf, 2, 1);
        end
        function e4e = get.e4e(mesh)
            nrF4e = size(mesh.fmask,2);
            e4e = zeros(nrF4e, mesh.nrE);
            e4e((mesh.e4fb-1)*nrF4e + mesh.ef4fb) = mesh.e4fb;
            e4e((mesh.e4f0-1)*nrF4e + mesh.ef4f0) = flipud(mesh.e4f0);
        end
        
        % Boundary facet affine information
        function J4fb  = get.J4fb(mesh),   if ~mesh.isFacetAnalyzed, mesh.affineFacets(); end, J4fb  = mesh.J4fb;  end
        function nx4fb = get.nx4fb(mesh),  if ~mesh.isFacetAnalyzed, mesh.affineFacets(); end, nx4fb = mesh.nx4fb; end
        function ny4fb = get.ny4fb(mesh),  if ~mesh.isFacetAnalyzed, mesh.affineFacets(); end, ny4fb = mesh.ny4fb; end
        function nz4fb = get.nz4fb(mesh),  if ~mesh.isFacetAnalyzed, mesh.affineFacets(); end, nz4fb = mesh.nz4fb; end
        
        % Interior facet affine information
        function J4f0  = get.J4f0(mesh),  if ~mesh.isFacetAnalyzed, mesh.affineFacets(); end, J4f0   = mesh.J4f0;  end
        function nx4f0 = get.nx4f0(mesh), if ~mesh.isFacetAnalyzed, mesh.affineFacets(); end, nx4f0  = mesh.nx4f0; end
        function ny4f0 = get.ny4f0(mesh), if ~mesh.isFacetAnalyzed, mesh.affineFacets(); end, ny4f0  = mesh.ny4f0; end
        function nz4f0 = get.nz4f0(mesh), if ~mesh.isFacetAnalyzed, mesh.affineFacets(); end, nz4f0  = mesh.nz4f0; end
        
        function tx4f0 = get.tx4f0(mesh), if mesh.dim ~= 2, error('Only 2D is supported'); end, if ~mesh.isFacetAnalyzed, mesh.affineFacets(); end, tx4f0  = -mesh.ny4f0; end
        function ty4f0 = get.ty4f0(mesh), if mesh.dim ~= 2, error('Only 2D is supported'); end, if ~mesh.isFacetAnalyzed, mesh.affineFacets(); end, ty4f0  = mesh.nx4f0; end
        
        function tx4fb = get.tx4fb(mesh), if mesh.dim ~= 2, error('Only 2D is supported'); end, if ~mesh.isFacetAnalyzed, mesh.affineFacets(); end, tx4fb  = -mesh.ny4fb; end
        function ty4fb = get.ty4fb(mesh), if mesh.dim ~= 2, error('Only 2D is supported'); end, if ~mesh.isFacetAnalyzed, mesh.affineFacets(); end, ty4fb  = mesh.nx4fb; end
        
        % Interior facet affine information
        function J4f  = get.J4f(mesh),  J4f   = [mesh.J4fb,  mesh.J4f0];  end
        function nx4f = get.nx4f(mesh), nx4f  = [mesh.nx4fb, mesh.nx4f0]; end
        function ny4f = get.ny4f(mesh), ny4f  = [mesh.ny4fb, mesh.ny4f0]; end
        function nz4f = get.nz4f(mesh), nz4f  = [mesh.nz4fb, mesh.nz4f0]; end
        
        function e4BC = e4BC(mesh, BCName), e4BC = mesh.e4fb(mesh.BC(BCName)); end
        function ef4BC = ef4BC(mesh, BCName), ef4BC = mesh.ef4fb(mesh.BC(BCName)); end
        
        %% DEPENDENT GET METHODS
        
        % Triangle coordinates
        function x4e = get.x4e(mesh), x4e = reshape(mesh.x(mesh.v4e), mesh.dim+1, []); end
        function y4e = get.y4e(mesh), y4e = reshape(mesh.y(mesh.v4e), mesh.dim+1, []); end
        function z4e = get.z4e(mesh), z4e = reshape(mesh.z(mesh.v4e), mesh.dim+1, []); end
        
        % Boundary facet coordinates
        function x4fb = get.x4fb(mesh), x4fb = reshape(mesh.x(mesh.v4fb), mesh.dim, []); end
        function y4fb = get.y4fb(mesh), y4fb = reshape(mesh.y(mesh.v4fb), mesh.dim, []); end
        function z4fb = get.z4fb(mesh), z4fb = reshape(mesh.z(mesh.v4fb), mesh.dim, []); end
        
        % Interior facet coordinates
        function x4f0 = get.x4f0(mesh), x4f0 = reshape(mesh.x(mesh.v4f0), mesh.dim, []); end
        function y4f0 = get.y4f0(mesh), y4f0 = reshape(mesh.y(mesh.v4f0), mesh.dim, []); end
        function z4f0 = get.z4f0(mesh), z4f0 = reshape(mesh.z(mesh.v4f0), mesh.dim, []); end
        
        % Interior facet coordinates
        function x4f = get.x4f(mesh), x4f = reshape(mesh.x(mesh.v4f), mesh.dim, []); end
        function y4f = get.y4f(mesh), y4f = reshape(mesh.y(mesh.v4f), mesh.dim, []); end
        function z4f = get.z4f(mesh), z4f = reshape(mesh.z(mesh.v4f), mesh.dim, []); end
        
        function x4BC = x4BC(mesh, name), x4BC = reshape(mesh.x(mesh.v4BC(name)), mesh.dim, []); end
        function y4BC = y4BC(mesh, name), y4BC = reshape(mesh.y(mesh.v4BC(name)), mesh.dim, []); end
        function z4BC = z4BC(mesh, name), z4BC = reshape(mesh.z(mesh.v4BC(name)), mesh.dim, []); end
        
        function v4BC  =  v4BC(mesh, name),  v4BC = mesh.v4fb(:,mesh.BC(name)); end
        function nx4BC = nx4BC(mesh, name), nx4BC = mesh.nx4fb (mesh.BC(name)); end
        function ny4BC = ny4BC(mesh, name), ny4BC = mesh.ny4fb (mesh.BC(name)); end
        function nz4BC = nz4BC(mesh, name), nz4BC = mesh.nz4fb (mesh.BC(name)); end
        function J4BC  =  J4BC(mesh, name),  J4BC = mesh.J4fb  (mesh.BC(name)); end
        function tx4BC = tx4BC(mesh, name), tx4BC = mesh.tx4fb (mesh.BC(name)); end
        function ty4BC = ty4BC(mesh, name), ty4BC = mesh.ty4fb (mesh.BC(name)); end
        
    end
    
    methods (Access = protected)
        
        affineSmplx(mesh);
        affineFacets(mesh);
        buildFacets(mesh);
        
    end
    
    methods (Access = protected)
        function propgrp = getPropertyGroups(mesh)
            if ~isscalar(mesh)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(fieldnames(mesh));
            else
                propList = cell(1, 2 + mesh.dim);
                propList(1:2) = {'shape', 'v4e'};
                d = 'xyz';
                for i = 1 : mesh.dim
                    propList{2 + i} = d(i);
                end
                
                propgrp = matlab.mixin.util.PropertyGroup(propList);
            end
        end
    end
end