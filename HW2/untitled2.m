function [x, xhist] = muller(f, x0, x1, x2, tol, maxit)

    if nargin < 6, maxit = 100; end
    if nargin < 5, tol   = 1e-12; end

    xhist = zeros(maxit+2,1);  
    xhist(1:3) = [x0; x1; x2];

    for k = 1:maxit
        f0 = f(x0);  f1 = f(x1);  f2 = f(x2);

        h0 = x1 - x0;
        h1 = x2 - x1;

        d0 = (f1 - f0) / h0;
        d1 = (f2 - f1) / h1;

        a = (d1 - d0) / (h1 + h0);
        b = a*h1 + d1;
        c = f2;

        disc = sqrt(b.^2 - 4*a*c);
        if abs(b + disc) > abs(b - disc)
            denom = b + disc;
        else
            denom = b - disc;
        end

        if abs(denom) == 0
            x3 = x2; 
        else
            x3 = x2 + (-2*c) / denom;
        end

        xhist(k+3) = x3;

        if abs(x3 - x2) <= tol*(1 + abs(x3)) || abs(f(x3)) <= tol
            xhist = xhist(1:k+3);
            x = x3;
            return;
        end

        x0 = x1;  x1 = x2;  x2 = x3;
    end

    x = x2;
    xhist = xhist(1:maxit+2);
end


f = @(x) x.^3 + x + 1;
[root1,it1] = muller(f,-1.0,0.0,1.0);

p = [1 0 1 1];
[q,~] = deconv(p,[1 -root1]);

a = q(1); b = q(2); c = q(3);
disc = b^2 - 4*a*c;
root2 = (-b + sqrt(disc))/(2*a);
root3 = (-b - sqrt(disc))/(2*a);
fprintf('root1 ≈ %.5e%+.5ei  (iters=%d)\n', real(root1), imag(root1), it1);
fprintf('root2 ≈ %.5e%+.5ei\n', real(root2), imag(root2));
fprintf('root3 ≈ %.5e%+.5ei\n', real(root3), imag(root3));
