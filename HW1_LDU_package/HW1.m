%% HW1 - LDU Decomposition (script version)
% This script demonstrates calling the LDU_Decomposition function.
% If you prefer a Live Script, open this file in MATLAB and choose:
%  Home -> Save As -> Live Script (.mlx)

% Example matrix
A = [2 1 1; 4 -6 0; -2 7 2];

% Call the function (make sure LDU_Decomposition.m is in the same folder)
[P, L, D, U] = LDU_Decomposition(A);

disp('P * A =');
disp(P*A);
disp('L * D * U =');
disp(L*D*U);

% Numerical verification
err = norm(P*A - L*D*U);
fprintf('Verification error = %.2e\n', err);
 