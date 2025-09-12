function [A, t] = transform(P, Q)

    p_avg = mean(P, 2);
    q_avg = mean(Q, 2);

    Pp = P - p_avg;
    Qp = Q - q_avg;

    B = Qp * Pp.';            

    [U, ~, V] = SVD(B);


    C = diag([1, 1, sign(det(U*V'))]);
    A = U * C * V.';

    t = q_avg - A * p_avg;
end
