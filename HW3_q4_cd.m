function HW3_q4_cd()
    clc; close all;

    % ===== Part (c): dominant plane via RANSAC =====
    file_c = 'cluttered_table.txt';   % adjust as needed
    fprintf('==== Part (c): %s ====\n', file_c);
    [pl_c, inliers_c, stats_c] = q4_c(file_c, 0.01, 4000);
    fprintf('RANSAC plane (normalized): a=%.6f, b=%.6f, c=%.6f, d=%.6f\n', ...
        pl_c.a, pl_c.b, pl_c.c, pl_c.d);
    fprintf('Inliers: %d / %d (%.1f%%)\n', stats_c.inliers, stats_c.N, 100*stats_c.inlier_ratio);
    fprintf('Mean |distance| on inliers: %.6f m\n\n', stats_c.avg_abs_dist);

    % ===== Part (d): four dominant planes =====
    file_d = 'clean_hallway.txt';     % adjust as needed
    fprintf('==== Part (d): %s ====\n', file_d);
    [planes_d, masks_d] = q4_d(file_d, 0.01, 4000, 200);
    fprintf('Found %d plane(s)\n', numel(planes_d));
    for i = 1:numel(planes_d)
        pl = planes_d{i};
        fprintf('  Plane %d: a=%.6f, b=%.6f, c=%.6f, d=%.6f\n', i, pl.a, pl.b, pl.c, pl.d);
    end
end

% ============================================================
function [plane, inlier_mask, stats] = q4_c(data_file, thresh, iters)
% Part (c): RANSAC dominant plane on a cluttered table
% Returns:
%   plane: struct with fields a,b,c,d (||[a b c]|| = 1, d signed)
%   inlier_mask: Nx1 logical for all points in data file
%   stats: struct with counts and mean inlier distance

    P = load_xyz(data_file);
    [plane, inlier_mask] = ransac_plane(P, thresh, iters);
    Pin = P(inlier_mask, :);
    avg_abs = mean(point_plane_abs_distance(Pin, plane));

    stats = struct();
    stats.N = size(P,1);
    stats.inliers = nnz(inlier_mask);
    stats.inlier_ratio = stats.inliers / stats.N;
    stats.avg_abs_dist = avg_abs;

    % Plot
    ttl = sprintf('Part (c): %s — dominant plane via RANSAC', data_file);
    plot_with_planes(P, {plane}, {inlier_mask}, ttl);
end

% ============================================================
function [planes, masks] = q4_d(data_file, thresh, iters, minCluster)
% Part (d): Extract four dominant planes iteratively (peel-off)
% Returns:
%   planes: 1xK cell of plane structs
%   masks : 1xK cell of Nx1 logical masks for original indexing

    if nargin < 4 || isempty(minCluster), minCluster = 200; end
    P = load_xyz(data_file);
    [planes, masks] = extract_k_planes(P, 4, thresh, iters, minCluster);
    ttl = sprintf('Part (d): %s — four dominant planes', data_file);
    plot_with_planes(P, planes, masks, ttl);
end

% ======================= RANSAC Core ========================
function [best_plane, best_mask] = ransac_plane(P, thresh, iters)
% RANSAC to find a dominant plane; refine with LS on inliers

    N = size(P,1);
    best_cnt = -1;
    best_plane = [];
    best_mask = false(N,1);

    for t = 1:iters
        idx = randperm(N,3);
        pl_cand = plane_from_3(P(idx(1),:), P(idx(2),:), P(idx(3),:));
        if isempty(pl_cand), continue; end

        d = point_plane_abs_distance(P, pl_cand);
        mask = d <= thresh;
        cnt  = nnz(mask);
        if cnt > best_cnt
            best_cnt   = cnt;
            best_plane = pl_cand;
            best_mask  = mask;
        end
    end

    if isempty(best_plane), error('RANSAC failed to find a plane.'); end

    % Refine by LS on inliers, keeping the representation consistent
    Pin = P(best_mask, :);
    pl_ref = fit_plane_ls_fixed_d(Pin, -1); % z = a x + b y + c => ax + by - z + c = 0
    best_plane = pl_ref;
end

function [planes, masks] = extract_k_planes(P, k, thresh, iters, minCluster)
% Iteratively peel off k dominant planes

    Rem = P;
    idx_map = (1:size(P,1)).';
    planes = {};
    masks  = {};

    for i = 1:k
        if size(Rem,1) < 3, break; end
        try
            [pl, mask_rem] = ransac_plane(Rem, thresh, iters);
        catch
            break
        end
        if nnz(mask_rem) < minCluster
            break
        end
        planes{end+1} = pl; %#ok<AGROW>

        % Map inliers back to original indexing
        in_global = false(size(P,1),1);
        in_global(idx_map(mask_rem)) = true;
        masks{end+1} = in_global; %#ok<AGROW>

        % Peel off
        Rem     = Rem(~mask_rem, :);
        idx_map = idx_map(~mask_rem);
    end
end

% ======================= Geometry Utils =====================
function P = load_xyz(path)
    P = load(path);
    if size(P,2) ~= 3
        error('Expected an N x 3 text file with columns [x y z].');
    end
end

function pl = fit_plane_ls_fixed_d(P, d_fixed)
% Fit ax + by + cz + d = 0 with fixed d (use c = -1 to match z = a x + b y + c)
% Here we choose c = -1 implicitly by fitting z = a x + b y + c first

    % Solve z = a*x + b*y + c via LS
    x = P(:,1); y = P(:,2); z = P(:,3);
    A = [x, y, ones(size(x))];
    t = A \ z;                % t = [a; b; c]
    a = t(1); b = t(2); c0 = t(3);

    % Convert to ax + by + cz + d = 0 form with c = -1, d = c0
    pl = normalize_plane(struct('a', a, 'b', b, 'c', -1, 'd', c0));
end

function pl = plane_from_3(p1, p2, p3)
% Plane through three points; normalized so ||[a b c]|| = 1

    v1 = p2 - p1; v2 = p3 - p1;
    n = cross(v1, v2);
    nrm = norm(n);
    if nrm < 1e-12
        pl = [];
        return;
    end
    n = n / nrm;
    d = -dot(n, p1);
    pl = normalize_plane(struct('a', n(1), 'b', n(2), 'c', n(3), 'd', d));
end

function pl = normalize_plane(pl)
% Normalize so ||(a,b,c)|| = 1 and flip sign so d <= 0 for consistency

    nrm = norm([pl.a pl.b pl.c]);
    if nrm == 0, return; end
    pl.a = pl.a/nrm; pl.b = pl.b/nrm; pl.c = pl.c/nrm; pl.d = pl.d/nrm;
    if pl.d > 0
        pl.a = -pl.a; pl.b = -pl.b; pl.c = -pl.c; pl.d = -pl.d;
    end
end

function d = point_plane_abs_distance(P, pl)
% |ax + by + cz + d| / ||n||
    num = abs(P * [pl.a; pl.b; pl.c] + pl.d);
    den = sqrt(pl.a^2 + pl.b^2 + pl.c^2);
    d = num / den;
end

% ======================= Visualization ======================
function plot_with_planes(P, planes, masks, ttl)
% P: Nx3, planes: {plane1, plane2, ...}, masks: {mask1, ...} or {}
    if nargin < 3, masks = {}; end

    figure; hold on; grid on; axis equal; view(45,25);

    plotted = false(size(P,1),1);
    if ~isempty(masks)
        C = lines(max(1, numel(masks)));
        for i = 1:numel(masks)
            m = masks{i};
            scatter3(P(m,1), P(m,2), P(m,3), 8, C(i,:), 'filled');
            plotted = plotted | m;
        end
    end
    if any(~plotted)
        scatter3(P(~plotted,1), P(~plotted,2), P(~plotted,3), 6, [0.6 0.6 0.6], 'filled');
    end

    % Extents for plane patches
    x_min = min(P(:,1)); x_max = max(P(:,1));
    y_min = min(P(:,2)); y_max = max(P(:,2));
    z_min = min(P(:,3)); z_max = max(P(:,3));
    xs = linspace(x_min, x_max, 28);
    ys = linspace(y_min, y_max, 28);
    zs = linspace(z_min, z_max, 28);

    for i = 1:numel(planes)
        pl = planes{i};
        a = pl.a; b = pl.b; c = pl.c; d = pl.d;
        if abs(c) >= max(abs(a),abs(b))
            [X,Y] = meshgrid(xs, ys);
            Z = (-(a*X + b*Y + d)) / c;
        elseif abs(b) >= max(abs(a),abs(c))
            [X,Z] = meshgrid(xs, zs);
            Y = (-(a*X + c*Z + d)) / b;
        else
            [Y,Z] = meshgrid(ys, zs);
            X = (-(b*Y + c*Z + d)) / a;
        end
        surf(X, Y, Z, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
    end

    xlabel('x (right)'); ylabel('y (down)'); zlabel('z (forward)');
    title(ttl);
end
