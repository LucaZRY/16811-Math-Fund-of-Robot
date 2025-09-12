function [P,L,D,U] = LDU_Decomposition(A)

    [m, n] = size(A);
    if m ~= n
        error('A must be a square matrix');
    end

    P = eye(m);
    L = eye(m);
    U = A;

    for k = 1:m-1
        [~, idx] = max(abs(U(k:m, k)));
        pivot_row = k + idx - 1;

        if U(pivot_row, k) == 0
            error('Matrix is singular at column %d.', k);
        end

        if pivot_row ~= k
            [P, U, L] = row_swap(U, P, L, pivot_row, k);
        end

        for i = k+1:m
            L(i,k) = U(i,k) / U(k,k);
            U(i,:) = U(i,:) - L(i,k) * U(k,:);
        end
    end

    D = diag(diag(U));
    U = D \ U;
end

function [P_new, U_new, L_new] = row_swap(U, P, L, r2, r1)
    U_new = U;
    P_new = P;
    L_new = L;

    U_new([r1 r2], :) = U_new([r2 r1], :);
    P_new([r1 r2], :) = P_new([r2 r1], :);

    if r1 > 1
        L_new([r1 r2], 1:r1-1) = L_new([r2 r1], 1:r1-1);
    end
end
