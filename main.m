clc
clear
close all

%% Parameteres
faceopaque = false;
icahedron_plot = true;

fs = 20e3;
c = 330; % m/s speed of sound

L = 3;
nPoints = 10*4^L + 2;
test_points_NUM = 2;

N = 10;
test_vel = 0.5*(rand(test_points_NUM, 3) - 0.5);
test_point_initial = 2*(rand(test_points_NUM, 3)-0.5);
test_point = zeros(test_points_NUM, 3, N);
test_point(:, :, 1) = test_point_initial;
for i = 2:N
    test_point(:, :, i) = test_point(:, :, i-1) + test_vel;
end

%% microphones location
theta = linspace(0, 2*pi, 9).';
mics_locs = 0.254/2*[cos(theta), sin(theta), zeros(length(theta), 1)];
mics_locs = mics_locs(1:end-1, :);
mics_locs2 = 0.254/2*[cos(theta), zeros(length(theta), 1), sin(theta)];
mics_locs2 = mics_locs2(1:end-1, :);

mics_locs = [mics_locs; mics_locs2];

%% creating test points
[vMat, fMat] = spheretri(nPoints);

%% plot

if icahedron_plot == true
    x = vMat(:, 1);
    y = vMat(:, 2);
    z = vMat(:, 3);
    if faceopaque == true
        h=trisurf(fMat,x,y,z, 'FaceColor', 'white', 'EdgeColor', 0.2*[1 1 1], 'LineWidth', 1, 'FaceAlpha', 1);
        hold on; p = plot3(x,y,z, 'ko', 'MarkerFaceColor', 'black');
        plot3(mics_locs(:, 1), mics_locs(:, 2), mics_locs(:, 3), 'sg')
    else
        h=trisurf(fMat,x,y,z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1],'LineWidth', 1 );
        hold on; p = plot3(x,y,z, 'ko', 'MarkerFaceColor', 'black');
        plot3(mics_locs(:, 1), mics_locs(:, 2), mics_locs(:, 3), 'sg', 'MarkerFaceColor', 'green')
    end
    
    xlabel('x-axis (m)');
    ylabel('y-axis (m)');
    zlabel('z-axis (m)');
    title( ['test points']);
    
    axis equal;
    grid off;
    
end

%% finding tao(differential time delay) for test points

tao_test = zeros(size(mics_locs, 1), size(mics_locs, 1), size(vMat, 1));

for i = 1:size(mics_locs, 1)
    for j = 1:size(mics_locs, 1)
        tao_test(i, j, :) = (fs/c)*dot(repmat(mics_locs(i, :) - mics_locs(j, :), size(vMat, 1), 1) , vMat, 2);
    end
end

% first mic is reference
steering_test = zeros(size(mics_locs, 1), size(vMat, 1));
for i = 1:size(mics_locs, 1)
    steering_test(i, :) = (fs/c)*dot(repmat(mics_locs(1, :) - mics_locs(i, :), size(vMat, 1), 1) , vMat, 2);
end

%%
for iter = 1:N
    test_point_norm = test_point(:, :, iter);
    test_point_norm = test_point_norm / norm(test_point_norm);
    
    hold on;plot3(test_point(:, 1, iter), test_point(:, 2, iter), test_point(:, 3, iter), 'rp', 'MarkerFaceColor', 'red')
    
    test_steer = zeros(size(mics_locs, 1), size(test_point, 1));
    
    % test_point = vMat(15, :);
    % hold on;plot3(test_point(:, 1), test_point(:, 2), test_point(:, 3), 'rp', 'MarkerFaceColor', 'red')
    for j = 1:size(test_point, 1)
        for i = 1:size(mics_locs, 1)
            test_steer(i, j) = (fs/c)*dot(mics_locs(1, :) - mics_locs(i, :) , test_point_norm(j, :), 2);
        end
    end
    
    for j = 1:size(test_point, 1)
        corr_test = repmat(test_steer(:, j), 1, size(steering_test, 2)) - steering_test;
        corr_test = sum(corr_test.^2).^0.5;
        % corr_test = dot(repmat(test_steer, 1, size(steering_test, 2)) , steering_test, 1);
        [~, ind_max(j)] = min(corr_test);
    end
    
    estimated_locs = [x(ind_max),y(ind_max),z(ind_max)];
    
    plot3(estimated_locs(:, 1), estimated_locs(:, 2), estimated_locs(:, 3), 'bv', 'MarkerFaceColor', 'blue')
    
end
%% tracking

% parameters
delta_T = 0.5;
F = eye(6);
for i = 4:6
    F(i-3, i) = delta_T;
end

sigma2_Q = 0.5;
Q = zeros(6);
for i = 4:6
    Q(i, i) =  sigma2_Q;
end

H = zeros(3, 6);
for i = 1:3
    H(i, i) = 1;
end

sigma2_R = 0.5;
R = sigma2_R*eye(3);







