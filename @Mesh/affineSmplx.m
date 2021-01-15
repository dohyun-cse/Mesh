function affineSmplx(mesh)

mesh.isTriAnalyzed = true;

% preallocate
if mesh.gpuFlag
    dxdr = zeros(mesh.dim, mesh.dim, mesh.nrE, 'gpuArray');
else
    dxdr = zeros(mesh.dim, mesh.dim, mesh.nrE);
end

% get dx/dr
ax = 'xyz';
for i = 1 : mesh.dim
    dxdr(i,:,:) = (mesh.(ax(i))(mesh.v4e(2:end,:)) - mesh.(ax(i))(mesh.v4e(1,:)))/2;
end

% get Jacobian using determinant of dx/dr
mesh.J = mydet(dxdr);
% get dr/dx by inverting dx/dr. We can utilize determinant to obtain
% inverse.
drdx = myinv(dxdr, mesh.J);

% reshape and store data with familar names.
mesh.J = reshape(mesh.J, 1, mesh.nrE);
phy_dir = 'xyz';
ref_dir = 'rsm';
for i = 1 : mesh.dim
    for j = 1 : mesh.dim
        mesh.([ref_dir(j) phy_dir(i)]) = reshape(drdx(j,i,:), 1, mesh.nrE);
    end
end
phy_dir = 'xyz';
ref_dir = 'rsm';
for i = 1 : mesh.dim
    for j = 1 : mesh.dim
        mesh.([phy_dir(i) ref_dir(j)]) = reshape(dxdr(i,j,:), 1, mesh.nrE);
    end
end
end

%% Vectorized determinant for 1x1, 2x2, 3x3
function d = mydet(A)
switch size(A,1)
    case 1
        d = A;
    case 2
        d = A(1,1,:).*A(2,2,:) - A(2,1,:).*A(1,2,:);
    case 3
        d = A(1,1,:).*A(2,2,:).*A(3,3,:) ...
          + A(1,2,:).*A(2,3,:).*A(3,1,:) ...
          + A(1,3,:).*A(2,1,:).*A(3,2,:) ...
          - A(1,1,:).*A(2,3,:).*A(3,2,:) ...
          - A(1,2,:).*A(2,1,:).*A(3,3,:) ...
          - A(1,3,:).*A(2,2,:).*A(3,1,:);
end     
end

%% Vectorized inverse for 1x1, 2x2, 3x3
function invA = myinv(A, d)
switch size(A,1)
    case 1
        invA = 1./A;
    case 2
        invA = [A(2,2,:), -A(1,2,:);
               -A(2,1,:),  A(1,1,:)]./reshape(d,1,1,[]);
    case 3
        invA = [A(2,2,:).*A(3,3,:)-A(2,3,:).*A(3,2,:), A(1,3,:).*A(3,2,:)-A(1,2,:).*A(3,3,:), A(1,2,:).*A(2,3,:)-A(1,3,:).*A(2,2,:);
                A(2,3,:).*A(3,1,:)-A(2,1,:).*A(3,3,:), A(1,1,:).*A(3,3,:)-A(1,3,:).*A(3,1,:), A(1,3,:).*A(2,1,:)-A(1,1,:).*A(2,3,:);
                A(2,1,:).*A(3,2,:)-A(2,2,:).*A(3,1,:), A(1,2,:).*A(3,1,:)-A(1,1,:).*A(3,2,:), A(1,1,:).*A(2,2,:)-A(1,2,:).*A(2,1,:)]./reshape(d,1,1,[]);
end
end