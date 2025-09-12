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

