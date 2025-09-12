function [U, S, V] = SVD(A)

    [m, n] = size(A);
    S = zeros(m, n);

    [U, D] = eig(A * A');
    [~, index] = sort(diag(D), 'descend');
    U = U(:, index);
    D = D(index, index);


    [V, Dv] = eig(A' * A);
    [~, index] = sort(diag(Dv), 'descend');
    V = V(:, index);
    Dv = Dv(index, index);

    r = rank(D);

    S(1:r, 1:r) = sqrt(D(1:r, 1:r));

    V(:, 1:r) = A' * U(:, 1:r) / S(1:r, 1:r);
end
