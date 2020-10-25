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

test_point = rand(1, 3)-0.5;
test_point = test_point / norm(test_point);
%% microphones location
theta = linspace(0, 2*pi, 9).';
mics_locs = 0.254/2*[cos(theta), sin(theta), zeros(length(theta), 1)];
mics_locs = mics_locs(1:end-1, :);

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
    
    hold on;plot3(test_point(:, 1), test_point(:, 2), test_point(:, 3), 'rp', 'MarkerFaceColor', 'red') 
    
    axis equal;
    view(2)
    
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

%% test

test_steer = zeros(size(mics_locs, 1), 1);

% test_point = vMat(15, :);
% hold on;plot3(test_point(:, 1), test_point(:, 2), test_point(:, 3), 'rp', 'MarkerFaceColor', 'red') 
for i = 1:size(mics_locs, 1)
    test_steer(i) = (fs/c)*dot(mics_locs(1, :) - mics_locs(i, :) , test_point, 2);
end

corr_test = repmat(test_steer, 1, size(steering_test, 2)) - steering_test;
corr_test = sum(corr_test.^2).^0.5;
% corr_test = dot(repmat(test_steer, 1, size(steering_test, 2)) , steering_test, 1);
[~, ind_max] = min(corr_test);

plot3(x(ind_max),y(ind_max),z(ind_max), 'bv', 'MarkerFaceColor', 'blue');












