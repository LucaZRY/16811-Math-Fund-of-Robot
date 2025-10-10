function HW3_4ab()

    clc; close all;

    data_file1 = 'clear_table.txt';     
    fprintf('==== Part (a): %s ====\n', data_file1);
    [t1, rss1, avg1] = q4_a_b(data_file1, true);
    fprintf('Plane (a): z = %.6f x + %.6f y + %.6f\n', t1(1), t1(2), t1(3));
    fprintf('RSS^0.5 (z residual): %.6f\n', rss1);
    fprintf('Mean point-to-plane distance: %.6f m\n\n', avg1);

    data_file2 = 'cluttered_table.txt';
    fprintf('==== Part (b): %s ====\n', data_file2);
    [t2, rss2, avg2] = q4_a_b(data_file2, true);
    fprintf('Plane (b): z = %.6f x + %.6f y + %.6f\n', t2(1), t2(2), t2(3));
    fprintf('RSS^0.5 (z residual): %.6f\n', rss2);
    fprintf('Mean point-to-plane distance: %.6f m\n', avg2);
end


function [t, rss, avg_dist] = q4_a_b(data_file, do_plot)


    P = load(data_file);    
    if size(P,2) ~= 3
        error('Expected an N x 3 text file with columns [x y z].');
    end
    x = P(:,1); y = P(:,2); z = P(:,3);

    A = [x, y, ones(length(x),1)];
    t = A \ z;                         

    zhat = A * t;
    rss  = sqrt(sum((z - zhat).^2));  

    a = t(1); b = t(2); c = t(3);
    numer = abs(a.*x + b.*y - z + c);
    denom = sqrt(a^2 + b^2 + 1);
    avg_dist = mean(numer ./ denom);

    if do_plot
        plot_points_and_plane(x, y, z, t, sprintf('%s — plane fit', data_file));
    end
end


function plot_points_and_plane(x, y, z, t, ttl)

    a = t(1); b = t(2); c = t(3);

    pad = 0.05;
    xs = linspace(min(x)-pad, max(x)+pad, 40);
    ys = linspace(min(y)-pad, max(y)+pad, 40);
    [X, Y] = meshgrid(xs, ys);
    Z = a.*X + b.*Y + c;

    figure; hold on; grid on; axis vis3d;
    scatter3(x, y, z, 8, 'r', 'filled');
    surf(X, Y, Z, 'FaceAlpha', 0.35, 'EdgeColor', 'none');

    xlabel('x (right)'); ylabel('y (down)'); zlabel('z (forward)');
    title(sprintf('%s\nz = %.4f x + %.4f y + %.4f', ttl, a, b, c));
    view(45, 25);
end
