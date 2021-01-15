function affineFacets(mesh)

% if facet mapping is not built,
if ~mesh.isFacetBuilt
    mesh.buildFacets();
end

% update flag
mesh.isFacetAnalyzed = true;

if mesh.dim == 1
    mesh.nx4fb = sign(diff(mesh.x4fb));
    mesh.J4fb = ones(1,mesh.nrfb);
    mesh.nx4f0 = sign(diff(mesh.x4f0));
    mesh.J4f0 = ones(1,mesh.nrf0);
else
    % preallocate
    aff4f = zeros(mesh.dim, mesh.dim-1, mesh.nrf, 'like', mesh.x);
    n4f = zeros(mesh.dim, mesh.nrf, 'like', mesh.x);

    % obtain dx/dr on each facet
    phy_dir = 'xyz';
    for i = 1 : mesh.dim
        aff4f(i,:,:) = (mesh.([phy_dir(i) '4f'])(2:end,:) - mesh.([phy_dir(i) '4f'])(1,:))/2;
    end

    % The normal direction can be obtained using determinant of minor matrix
    targ = true(mesh.dim,1);
    for i = 1 : mesh.dim
        targ(:) = true; targ(i) = false;
        if mod(i, 2)
            n4f(i,:) = mydet(aff4f(targ,:,:));
        else
            n4f(i,:) = -mydet(aff4f(targ,:,:));
        end
    end

    % Jacobian is the size of normal vector
    J4f = vecnorm(n4f,2);
    mesh.J4fb = J4f(1:mesh.nrfb);
    mesh.J4f0 = J4f(mesh.nrfb + 1:mesh.nrf);

    % Normalize
    n4f = n4f./J4f;

    % stoare data in familiar names
    for i = 1 : mesh.dim
        mesh.(['n' phy_dir(i) '4fb']) = n4f(i,1:mesh.nrfb);
        mesh.(['n' phy_dir(i) '4f0']) = n4f(i,mesh.nrfb+1:end);
    end
end
end

%% Vectroized determinant for 1x1, 2x2, and 3x3.
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