function plot_points(P, Q, A, t)
% PLOT_POINTS  Visualize P (blue), Q (red), and A*P+t (green). Assumes 3x3.

    p_avg = mean(P,2);  q_avg = mean(Q,2);
    R = A*P + t;        r_center = A*p_avg + t;

    figure; hold on; grid on; axis equal;
    plot3(P(1,:),P(2,:),P(3,:),'b.','MarkerSize',18);
    plot3(Q(1,:),Q(2,:),Q(3,:),'r.','MarkerSize',18);
    plot3(R(1,:),R(2,:),R(3,:),'g.','MarkerSize',18);

    plot3(p_avg(1),p_avg(2),p_avg(3),'b*','MarkerSize',8);
    plot3(q_avg(1),q_avg(2),q_avg(3),'r*','MarkerSize',8);
    plot3(r_center(1),r_center(2),r_center(3),'go','MarkerSize',8);

    for k = 1:3
        plot3([Q(1,k) R(1,k)], [Q(2,k) R(2,k)], [Q(3,k) R(3,k)], 'k-');
    end

    legend('P','Q','A P + t','centroid P','centroid Q','centroid AP+t');
    xlabel X; ylabel Y; zlabel Z; title('Rigid alignment (3×3) via SVD');
    hold off;
end
