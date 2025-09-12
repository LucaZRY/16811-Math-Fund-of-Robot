

function [P,L,D,U] = LDU(A)
    [m, n] = size(A);

    if m ~= n
        error('A must be a square matrix')
    end


    L = eye(m);
    P = eye(m);
    A_copy = A;

    for i = 1:m
        if A_copy(i,i) == 0
            change_row = i;
            for j = i+1:m
                if A_copy(j,i) ~= 0
                    change_row = j;
                    break
                end
            end
        end
        if change_row ~= i
            [P, A_copy, L] = change_row(A_copy, P, L, change_row, i);
        else
            error('Matrix is singular')
        end

        for j = i+1:m
            factor = A_copy(j,i) / A_copy(i,i);
            A_copy(j,:) = A_copy(j,:) - factor * A_copy(i,:);
            L(j,i) = factor;
        end
    end

    D = diag(diag(A_copy));
    U = D \ A_copy;
end

function [P_new, A_new, L_new] = change_row(A, P, L, row2, row1)
%ROW_SWAP Swap two rows in A, P, and adjust L accordingly.
%
%   Inputs:
%       A    - Matrix to modify
%       P    - Permutation matrix
%       L    - Lower-triangular matrix
%       row2 - Row to swap into row1
%       row1 - Target row index
%
%   Outputs:
%       P_new - Updated permutation matrix
%       A_new - Updated A matrix
%       L_new - Updated lower-triangular matrix

    A_new = A;
    P_new = P;
    L_new = L;

    % Swap rows in A and P
    A_new([row1,row2], :) = A_new([row2,row1], :);
    P_new([row1,row2], :) = P_new([row2,row1], :);

    % Swap only the first (row1-1) columns in L
    if row1 > 1
        L_new([row1,row2], 1:row1-1) = L_new([row2,row1], 1:row1-1);
    end
end

