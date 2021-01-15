function buildFacets(mesh)

% update flag
mesh.isFacetBuilt = true;

% the number of vertices for one facet
nrV4f = size(mesh.fmask,1);
% the number of facets for one triangle
nrF4e = size(mesh.fmask,2);

% all facets surrounding each triangle, with duplication
v4F = reshape(mesh.v4e(mesh.fmask,:), nrV4f, []);
% remove duplication and store mapping between unique list and original
% list
[~, ind, bak] = unique(sort(v4F).', 'rows');
% restore original ordering so that the resulting facet vertices are 
% ordered in a consistent direction.
v4F = v4F(:,ind);

% initialize index as 2
ii = repmat(2, mesh.nrE*nrF4e, 1);
% If the facet is first one appeared in the unique list, mark as 1
% those will be the positive element of a facet.
ii(ind) = 1;

% the number of facets
nrF = size(v4F, 2);

% preallocate neighboring triangle information
e4F = zeros(2, nrF);
% Recall that original facet lists are built triangle-wise,
% the neighboring triangle's index sucessively increases
e4F(ii + (bak-1)*2) = reshape(repmat(1:mesh.nrE, nrF4e, 1), [], 1);

% Again, neighboring triangles' local facet index increases succesively 
ef4F = zeros(2, nrF);
ef4F(ii + (bak-1)*2) = repmat(1:nrF4e, 1, mesh.nrE);

% If there is no negative element, it is boundary
isBoundary = e4F(2,:) == 0;

% Store data
mesh.v4fb = v4F(:, isBoundary);
mesh.v4f0 = v4F(:,~isBoundary);

mesh.e4fb = e4F(1, isBoundary); mesh.ef4fb = ef4F(1, isBoundary);
mesh.e4f0 = e4F(:,~isBoundary); mesh.ef4f0 = ef4F(:,~isBoundary);
if mesh.dim == 3
    [i, ~] = find(mesh.v4f0(1,:) == mesh.v4e(mesh.fmask_neg(:,mesh.ef4f0(2,:)) + 4*(mesh.e4f0(2,:)-1)));
    mesh.ef4f0(2,:) = mesh.ef4f0(2,:) + (i.'-1)*4;
end