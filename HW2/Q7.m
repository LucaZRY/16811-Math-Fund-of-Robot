clear; clc;

syms x y real
p = 2*x^2 + 2*y^2 - 4*x - 4*y + 3;
q = x^2 + y^2 + 2*x*y - 5*x - 3*y + 4;

%% (a) Plot zero-contours p=0 and q=0
figure('Color','w'); hold on; grid on;
fimplicit(p==0, [-1 3 -1 3], 'LineWidth',1.8);
fimplicit(q==0, [-1 3 -1 3], 'LineWidth',1.8);
axis equal; xlabel('x'); ylabel('y');
legend('p(x,y)=0','q(x,y)=0','Location','best');
title('Zero-contours of p and q');