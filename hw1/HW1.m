A = [ 2.0  1.0  1.0;
       4.0 -6.0  0.0;
      -2.0  7.0  2.0 ];
[P, L, D, U] = LDU(A);

% Display results
disp('Permutation matrix P:');
disp(P);

disp('Lower triangular matrix L:');
disp(L);

disp('Diagonal matrix D:');
disp(D);

disp('Upper triangular matrix U:');
disp(U);

% Check if PA = LDU
fprintf('Error = %.2e\n', norm(P*A - L*D*U, 'fro'));