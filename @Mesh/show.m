function h = show(mesh, varargin)
% Mesh display routines
args = varargin;
if ~isempty(args) && isa(args{1}, 'matlab.graphics.axis.Axes')
    ax = args{1};
    args(1) = [];
else
    ax = gca;
end
if ~ishold(ax)
    cla(ax, 'reset');
end

switch mesh.dim
    case 1
        H = plot(ax, reshape([mesh.x4e;NaN(1, mesh.nrE)],[],1), zeros(3*mesh.nrE,1), '-kd', args{:});
    case 2
        H = patch(ax, 'Faces', mesh.v4e.', 'Vertices', [mesh.x; mesh.y].', 'FaceColor', [0.7, 0.7, 0.7], args{:});
        if nargout
            h = H;
        end
    case 3
        
        H = patch(ax, ...
            'XData', reshape(mesh.x4e(mesh.fmask,:), 3, []), ...
            'YData', reshape(mesh.y4e(mesh.fmask,:), 3, []), ...
            'ZData', reshape(mesh.z4e(mesh.fmask,:), 3, []), ...
            'FaceColor', [0.5, 0.5, 0.5], ...
            'FaceAlpha', 0.2, ...
            args{:});
end


if ~ishold(ax)
    set(ax, 'Color', 'none');
    axis(ax, 'tight', 'equal');
    if mesh.dim == 1
        set(ax,'YLim', ([1,-1;-1,1]*get(gca, 'XLim').'/2).');
        box(ax, 'off');
        set(ax, 'YTick', []);
    end
    if mesh.dim == 3
        view(ax, 3);
    end
else
    hold(ax, 'off');
end
if nargout
    h = H;
end
end