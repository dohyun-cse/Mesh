# Mesh (MATLAB class)
1. Supports 1D, 2D, 3D simplices (line, triangle, tetrahedron)
2. Mesh display function
3. Automatic facet generation and neighboring information
4. Automatic generation of Affine information

See `example.mlx` for more information.

__New issues and contributions are welcomed!__
## How to use?
You can create a `Mesh` object by
```MATLAB
    mesh = Mesh(v4e, x); % 1D line mesh / v4e: 2 x nrE
    mesh = Mesh(v4e, x, y); % 2D triangular mesh / v4e: 3 x nrE
    mesh = Mesh(v4e, x, y, z); % 3D tetrahedral mesh / v4e: 4 x nrE
```
where `v4e(:,i)` contains vertex index for each element.

Or, you can get uniform mesh by
```MATLAB
    mesh = line2linemesh(0,1,10); % 1D uniform mesh in (0,1) with 10 elements.
    mesh = rect2rectmesh(0,1,0,2,10,20); % 2D uniform mesh in (0,1)x(0,2) with 10 x 20 (x2) elements.
    mesh = cube2cubemesh(0,1,0,2,0,3,10,20,30); % 1D uniform mesh in (0,1)x(0,2)x(0,3) with 10 x 20 x 30 (x6) elements.
```

## Display Mesh
You can use `show` method to display mesh

```MATLAB
    mesh.show(); % default configuration
    mesh.show('FaceColor', 'r'); % change face color
    mesh.show(ax, __); % specify target axes
```
Here, the optional arguments are from `plot` for 1D and `patch` for 2D and 3D.



## Facet Information
When a `Mesh` is created, facet information can be obtained, where
- `v4f` is the __vertex__ index for each facet
- `x4f` is the __x-coordinate__ for each facet (or y, z)
- `nx4f` is the __facet normal__ for each facet (or y, z)
- `J4f` is the __Jacobian__ for each facet
- `e4f` is the __neighboring elements' index__ for each facet
- `ef4f` is the __neighboring elements' local facet index__ for each facet
- `fmask` and `fmask_neg` are the local facet vertex index of an element (positive / negative)

All information about facets can be obtained seperately for interior/boundary facets  
by changing `f` to `f0` or `f0`.

If you want to access to the first neighboring element to 10th facet, you can use
```MATLAB
    mesh.e4f(1,10); % gives the first neighboring element of 10th facet
    mesh.ef4f(1,10); % gives the positive element's which facet is 10th facet.
```
In other words,
```MATLAB
    mesh.v4f(:,10) == mesh.v4e(mesh.fmask(:,mesh.ef4f(1,10)), mesh.e4f(1,10))
```
Or, for negative element (the second neighboring element)
```MATLAB
    mesh.v4f(:,10) == mesh.v4e(mesh.fmask_neg(:,mesh.ef4f(2,10)), mesh.e4f(2,10))
```

## Affine Information
In finite element assembly, affine information of physical triangles is often used.

One can access to the affine information by
```MATLAB
    % dX/dR
    mesh.xr, mesh.xs, mesh.xm
    mesh.yr, mesh.ys, mesh.ym
    mesh.zr, mesh.zs, mesh.zm

    % Jacobial
    mesh.J

    % dR/dX
    mesh.rx, mesh.ry, mesh.rz
    mesh.sx, mesh.sy, mesh.sz
    mesh.mx, mesh.my, mesh.mz
```