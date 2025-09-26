f  = @(x) tan(x) - x;     
fp = @(x) tan(x).^2;      

x = linspace(13, 17.5, 1000);
y = f(x);

figure;
plot(x, y, 'b-', 'LineWidth', 1.5); hold on;
yline(0, 'k--');     
xline(15, 'r--');    
ylim([-20 20]);      
xlabel('x');
ylabel('f(x) = tan(x) - x');
title('f(x) = tan(x) - x with roots near x=15');

x_low  = newton(f, fp, 14);
x_high = newton(f, fp, 17.25);

fprintf('x_low  = %.15f\n', x_low);
fprintf('x_high = %.15f\n', x_high);

function root = newton(f, fp, x0, tol, max_iter)
    if nargin < 4, tol = 1e-12; end
    if nargin < 5, max_iter = 50; end
    
    x = x0;
    for k = 1:max_iter
        x_new = x - f(x)/fp(x);
        if abs(x_new - x) < tol
            root = x_new;
            return;
        end
        x = x_new;
    end
    error('Newton did not converge');
end