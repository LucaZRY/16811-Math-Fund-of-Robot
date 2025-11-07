function hw3_q4_alt()
clc; close all;

data = load('clear_table.txt');  
x = data(:,1); y = data(:,2); z = data(:,3);
P = [x y z];
planeA = fit_plane_pca(P); % get [a b c d]
planeA = normalize_plane(planeA);  
[mean_absA, rmsA] = plane_stats(P, planeA);

fprintf('4.a PCA plane (ax+by+cz+d=0): [%.6g  %.6g  %.6g  %.6g]\n', planeA);
fprintf('       mean|dist|=%.6g,  RMS=%.6g\n\n', mean_absA, rmsA);

figure('Name','4.a PCA on clear_table'); hold on; grid on; axis equal;
scatter3(x, z, y, 8, 'b', 'filled');
draw_plane_xz(planeA, x, z, y, 'FaceAlpha',0.35, 'FaceColor',[0.9 0.4 0.1], 'EdgeColor','none');
xlabel('x [m]'); ylabel('z [m]'); zlabel('y [m]');
title('4.a PCA Fit (x–z base, y vertical)'); view(3);

%% -------------------------- 4.b (SVD TLS) ----------------------
data = load('cluttered_table.txt');
x = data(:,1); y = data(:,2); z = data(:,3);
P = [x y z];

% Total least squares by SVD on [x y z 1] (NO robustness) — expect bias on clutter.
planeB = fit_plane_tls_homog(P);
planeB = normalize_plane(planeB);
[mean_absB, rmsB] = plane_stats(P, planeB);
fprintf('4.b TLS (SVD) plane on clutter: [%.6g  %.6g  %.6g  %.6g]\n', planeB);
fprintf('       mean|dist|=%.6g,  RMS=%.6g (note: biased if strong outliers)\n\n', mean_absB, rmsB);

figure('Name','4.b TLS (SVD) on cluttered_table'); hold on; grid on; axis equal;
scatter3(x, z, y, 8, [0.6 0.6 0.6], 'filled');
draw_plane_xz(planeB, x, z, y, 'FaceAlpha',0.3, 'FaceColor',[0.2 0.6 1.0], 'EdgeColor','none');
xlabel('x [m]'); ylabel('z [m]'); zlabel('y [m]');
title('4.b TLS (no robustness) — visibly biased by outliers'); view(3);

%% -------------------------- 4.c (IRLS) -------------------------
% Robust single-plane fit via IRLS (Huber). This replaces your RANSAC.
opts.delta = 0.02;      % Huber delta in meters (tune)
opts.max_iters = 50;
opts.verbose = true;
[planeC, wC] = robust_plane_irls(P, opts);
planeC = normalize_plane(planeC);
[mean_absC, rmsC] = plane_stats(P, planeC);
fprintf('4.c IRLS (Huber) plane on clutter: [%.6g  %.6g  %.6g  %.6g]\n', planeC);
fprintf('       mean|dist|=%.6g,  RMS=%.6g\n\n', mean_absC, rmsC);

figure('Name','4.c IRLS on cluttered_table'); hold on; grid on; axis equal;
% visualize weights (higher weight = more likely inlier)
scatter3(x, z, y, 10, wC, 'filled'); colormap turbo; colorbar; caxis([0 1]);
draw_plane_xz(planeC, x, z, y, 'FaceAlpha',0.25, 'FaceColor',[0.1 0.8 0.3], 'EdgeColor','none');
xlabel('x [m]'); ylabel('z [m]'); zlabel('y [m]');
title('4.c Robust IRLS (weights show inlier likelihood)'); view(3);

%% -------------------------- 4.d (multi-plane, clean) -----------
data = load('clean_hallway.txt');
x = data(:,1); y = data(:,2); z = data(:,3);
P = [x y z];

K = 4;                        % find up to 4 planes
opts.delta = 0.01;            % tighter delta in clean data
opts.max_iters = 40; opts.verbose = false;

remain = true(size(P,1),1);
planesD = [];
colors = lines(K);
figure('Name','4.d Multi-plane (clean hallway)'); hold on; grid on; axis equal;
scatter3(x, z, y, 3, [0.7 0.7 0.7], 'filled');

for k = 1:K
    if nnz(remain) < 3, break; end
    [plane, w] = robust_plane_irls(P(remain,:), opts);
    plane = normalize_plane(plane);
    planesD = [planesD; plane];

    % inliers by robust sigma (MAD) around this plane
    d = point_to_plane_signed(P(remain,:), plane);
    sig = 1.4826 * mad(d, 1);
    thr = max(2.5*sig, 0.01);
    inliers_local = abs(d) <= thr;

    % report stats
    [mn, rm] = plane_stats(P(remain,:), plane);
    fprintf('4.d Plane %d: [%.3e %.3e %.3e %.3e] | mean|dist|=%.3e | RMS=%.3e | inliers=%d\n', ...
        k, plane, mn, rm, nnz(inliers_local));

    % plot this plane
    draw_plane_xz(plane, x, z, y, 'FaceAlpha',0.3, 'FaceColor',colors(k,:), 'EdgeColor','none');

    % peel off inliers
    idx = find(remain);
    remain(idx(inliers_local)) = false;
end
xlabel('x [m]'); ylabel('z [m]'); zlabel('y [m]');
title('4.d Clean Hallway: Sequential IRLS + peel-off'); view(3);

%% -------------------------- 4.e (multi-plane, clutter) ---------
data = load('cluttered_hallway.txt');
x = data(:,1); y = data(:,2); z = data(:,3);
P = [x y z];

K = 4;                         % up to 4 dominant planes
opts.delta = 0.02;             % slightly looser for clutter
opts.max_iters = 60; opts.verbose = false;

remain = true(size(P,1),1);
planesE = [];
figure('Name','4.e Multi-plane (cluttered hallway)'); hold on; grid on; axis equal;
scatter3(x, z, y, 3, [0.7 0.7 0.7], 'filled');

avg_smooth = [];
for k = 1:K
    if nnz(remain) < 3, break; end
    [plane, w] = robust_plane_irls(P(remain,:), opts);
    plane = normalize_plane(plane);

    d_all = point_to_plane_signed(P(remain,:), plane);
    sig = 1.4826 * mad(d_all, 1);
    thr_in    = max(2.0*sig, 0.01);  % inlier acceptance
    thr_peel  = max(3.0*sig, 0.02);  % peel a bit wider
    inliers_local = abs(d_all) <= thr_in;

    [mn, rm] = plane_stats(P(remain,:), plane);
    planesE = [planesE; plane];
    avg_smooth(k,1) = mn;

    fprintf('4.e Plane %d: [%.3e %.3e %.3e %.3e] | mean|dist|=%.3e | RMS=%.3e | inliers=%d\n', ...
        k, plane, mn, rm, nnz(inliers_local));

    draw_plane_xz(plane, x, z, y, 'FaceAlpha',0.28, 'EdgeColor','none');

    % peel-off using the larger threshold (keeps dominance, increases robustness)
    idx = find(remain);
    peel_mask = abs(d_all) <= thr_peel;
    remain(idx(peel_mask)) = false;
end

% Report the "smoothest" (smallest mean|dist|)
if ~isempty(avg_smooth)
    [best_val, safest_idx] = min(avg_smooth);
    fprintf('\n4.e Smoothest surface: Plane %d  |  mean|dist| = %.4e m\n', safest_idx, best_val);
end

xlabel('x [m]'); ylabel('z [m]'); zlabel('y [m]');
title('4.e Cluttered Hallway: Robust IRLS + MAD peel-off'); view(3);

end % <-- end main


%% ==================== Helper Functions ====================
function plane = fit_plane_pca(P)

    mu = mean(P, 1);
    Q  = P - mu;
    C  = (Q.' * Q) / size(Q,1); 
    [V, D] = eig(C, 'vector');
    [~, idx] = min(D);
    n = V(:, idx);
    d = -dot(n, mu.');
    plane = [n(:).' d];
end

function plane = fit_plane_tls_homog(P)
% Homogeneous TLS plane: smallest right singular vector of [x y z 1].
    X = [P, ones(size(P,1),1)];
    [~,~,V] = svd(X, 'econ');
    v = V(:, end);
    plane = v.';   % [a b c d]
end

function plane = normalize_plane(plane)
% Normalize to ||[a b c]|| = 1 and fix a sign convention (c<=0).
    n = plane(1:3);
    s = norm(n);
    if s < eps, error('Degenerate plane normal.'); end
    plane = plane / s;
    if plane(3) > 0, plane = -plane; end
end

function d = point_to_plane_signed(P, plane)
% Signed distances given normalized plane [a b c d].
    a=plane(1); b=plane(2); c=plane(3); dd=plane(4);
    d = a*P(:,1) + b*P(:,2) + c*P(:,3) + dd;
end

function [mean_abs, rmsd] = plane_stats(P, plane)
% Mean absolute and RMS signed distances.
    d = point_to_plane_signed(P, plane);
    mean_abs = mean(abs(d));
    rmsd = sqrt(mean(d.^2));
end

function [plane, w] = robust_plane_irls(P, opts)
% IRLS with Huber loss for plane fitting (total least squares style):
%  - Start with PCA/TLS to get an initial plane
%  - Alternate: compute residuals -> Huber weights -> weighted covariance -> update normal
%  - Return plane normalized, and final weights (0..1)
%
% opts.delta: Huber delta (meters)
% opts.max_iters
% opts.verbose

    delta = getfielddef(opts, 'delta', 0.02);
    max_iters = getfielddef(opts, 'max_iters', 50);
    verbose   = getfielddef(opts, 'verbose', false);

    % init with PCA plane
    plane = fit_plane_pca(P);
    plane = normalize_plane(plane);

    N = size(P,1);
    w = ones(N,1);

    for it = 1:max_iters
        d = point_to_plane_signed(P, plane);
        w_prev = w;

        % Huber weights based on residuals
        absd = abs(d);
        w = ones(N,1);
        idx = absd > delta;
        w(idx) = delta ./ absd(idx);

        % Weighted centroid and covariance
        W = w / sum(w);
        mu = sum(P .* W, 1);
        Q = P - mu;
        % Weighted covariance (diagonal weights)
        C = (Q.' * (Q .* W)) ;  % sum_i w_i (q_i q_i^T)

        % Update normal as smallest eigenvector
        [V,D] = eig(C, 'vector');
        [~, j] = min(D);
        n = V(:, j);
        d0 = -dot(n, mu.');

        plane = normalize_plane([n(:).' d0]);

        if verbose && mod(it,10)==0
            fprintf('    IRLS iter %d: mean|d|=%.4g, median|d|=%.4g\n', it, mean(absd), median(absd));
        end

        if max(abs(w - w_prev)) < 1e-4
            break;
        end
    end
end

function val = getfielddef(s, name, def)
    if isfield(s, name), val = s.(name); else, val = def; end
end

function draw_plane_xz(plane, x, z, y, varargin)
    % Robust plane drawer that avoids dividing by tiny coefficients.
    % Keeps x–z on the base axis when possible, but switches bases if needed.

    a=plane(1); b=plane(2); c=plane(3); d=plane(4);

    xr = linspace(min(x), max(x), 30);
    zr = linspace(min(z), max(z), 30);
    yr = linspace(min(y), max(y), 30);

    % choose which variable to solve for: the one with the LARGEST |coef|
    [~, idx] = max(abs([a b c]));  % 1->x, 2->y, 3->z

    switch idx
        case 2  % |b| largest -> solve for y (preferred x–z base)
            [X,Z] = meshgrid(xr, zr);
            Y = -(a*X + c*Z + d)/b;
            % clamp extreme Y to keep picture readable
            lo = prctile(y, 2); hi = prctile(y, 98);
            Y = max(min(Y, hi), lo);
            surf(X, Z, Y, varargin{:});

        case 3  % |c| largest -> solve for z (switch to x–y base)
            [X,Y] = meshgrid(xr, yr);
            Z = -(a*X + b*Y + d)/c;
            % keep axes order consistent: x (X), z (Z), y (Y)
            surf(X, Z, Y, varargin{:});

        case 1  % |a| largest -> solve for x (switch to y–z base)
            [Y,Z] = meshgrid(yr, zr);
            X = -(b*Y + c*Z + d)/a;
            % keep axes order consistent: x (X), z (Z), y (Y)
            surf(X, Z, Y, varargin{:});
    end
end