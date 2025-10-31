%% create obstacle field
close all
clear all
waypoints=300;
N=101;
OBST = [20,30;60,40;70,85];
epsilon = [25; 20; 30];

obs_cost = double(zeros(N));
for i=1:size(OBST,1)

    t = zeros(N);
    t(OBST(i,1),OBST(i,2)) = 1; %point obstacles
    
    t_cost = double(bwdist(t));
    t_cost(t_cost>epsilon(i))=epsilon(i);
    t_cost = 1/(2*epsilon(i))*(t_cost-epsilon(i)).^2;
    
    obs_cost = obs_cost + t_cost(1:N, 1:N);
end

figure(1)
imagesc(obs_cost')
% obstacle cost gradient
gx = diff(double(obs_cost),1,1);
gy = diff(double(obs_cost),1,2);
hold on


figure(1);
surface(1:N,1:N,double(obs_cost'));
xlabel('X')
hold on;

%% initial path
%world params
SX = 10; % START
SY = 10;
GX = 90; % GOAL
GY = 90;


traj = zeros(2,waypoints);
traj(1,1)=SX;
traj(2,1)=SY;
dist_x = GX-SX;
dist_y = GY-SY;
for i=2:waypoints
    traj(1,i)=traj(1,i-1)+dist_x/(waypoints-1);
    traj(2,i)=traj(2,i-1)+dist_y/(waypoints-1);
end

path_init = traj';
tt=size(path_init,1);
path_init_values = zeros(size(path_init,1),1);
for i=1:tt
    path_init_values(i)=obs_cost(floor(path_init(i,1)),floor(path_init(i,2)));
end
plot3(path_init(:,1),path_init(:,2),path_init_values,'.r','MarkerSize',20);
hold on

path = path_init;


%% Optimize it...
% your code comes here

path_b = path_init;

w_obs = 0.8;    % weight for obstacle cost gradient
w_smo = 4.0;    % weight for smoothness (xi - xi-1)
alpha = 0.1;    % step-size scaling
max_iter = 500;

for iter = 1:max_iter

    for i = 2:tt-1
        x = path_b(i,1); y = path_b(i,2);

        if x < 1 || x > N-1 || y < 1 || y > N-1
            continue;
        end

        ix = fix(x); iy = fix(y);
        interp_x1 = (x - ix) * (gx(ix+1, iy)   - gx(ix, iy))   + gx(ix, iy);
        interp_x2 = (x - ix) * (gx(ix+1, iy+1) - gx(ix, iy+1)) + gx(ix, iy+1);
        grad_x = (y - iy) * (interp_x2 - interp_x1) + interp_x1;

        interp_y1 = (x - ix) * (gy(ix+1, iy)   - gy(ix, iy))   + gy(ix, iy);
        interp_y2 = (x - ix) * (gy(ix+1, iy+1) - gy(ix, iy+1)) + gy(ix, iy+1);
        grad_y = (y - iy) * (interp_y2 - interp_y1) + interp_y1;

        grad_obs = [grad_x, grad_y];  

        grad_smooth = (path_b(i,:) - path_b(i-1,:)); 

        total_step = - w_obs * grad_obs - w_smo * grad_smooth;

        path_b(i,:) = path_b(i,:) + alpha * total_step;
    end

    % quick snapshots
    if iter == 100 || iter == 200 || iter == 500
        pv = zeros(tt,1);
        for k = 1:tt
            pv(k) = obs_cost(floor(path_b(k,1)), floor(path_b(k,2)));
        end
        figure(10); if iter==100, clf; end
        surf(1:N,1:N,double(obs_cost')), view(-25,70); hold on;
        xlabel('X'); xlim([0 100]); ylim([0 100]);
        plot3(path_init(:,1), path_init(:,2), path_init_values, '.r', 'MarkerSize', 12);
        plot3(path_b(:,1),   path_b(:,2),   pv,                 '.g', 'MarkerSize', 18);
        title(sprintf('Part (b): one-sided smoothness, iter = %d', iter));
        hold off;
        drawnow;
    end
end



%% plot the trajectories
path_values = zeros(tt,1);
for i=1:tt
    path_values(i)=obs_cost(floor(path(i,1)),floor(path(i,2)));
end
figure(1)
hold on;
plot3(path(:,1),path(:,2),path_values,'.g','MarkerSize',20);

hold off;
