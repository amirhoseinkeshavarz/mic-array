% clc
clear
close all

%% Parameteres
faceopaque = false;
icahedron_plot = true;

[voice1 , fs] = audioread('starter.wav');
[voice2 , ~] = audioread('heart.wav');
[voice3 , ~] = audioread('footsteps.wav');
time = 10; % second
nSample = time * fs;
voice48K = [voice1(1:nSample) voice2(1:nSample) voice3(1:nSample)];
ds = 3;
voice16K = voice48K(1:ds:end,:); % downsample from 48KS/s to 16 KS/s
fs = fs/ds;
% sound(voice16K(:,1),fs)

c = 343; % m/s speed of sound

L = 3;
spacePoints = 10*4^L + 2;

%% microphones location
theta = linspace(0, 2*pi, 9).';
mics_locs = 0.254/2*[cos(theta), sin(theta), zeros(length(theta), 1)];
mics_locs = mics_locs(1:end-1, :);
% mics_locs2 = 0.254/2*[cos(theta), zeros(length(theta), 1), sin(theta)];
% mics_locs2 = mics_locs2(1:end-1, :);

mics_locs = [mics_locs];

%% creating test points
[vMat, fMat] = spheretri(spacePoints);
aa = vMat(:,3)<0;
vMat(aa,:) = [];

%% Target Modelling
N = 1;
% test_vel = [0.5*(rand(test_points_NUM, 2) - 0.5) rand(test_points_NUM, 1)];
targetIndx = [120  ];
P = [1 1 1];
test_point = vMat(targetIndx,:);

Display

%% finding tau(differential time delay) for test points

tau_test = zeros(size(mics_locs, 1), size(mics_locs, 1), size(vMat, 1));

for i = 1:size(mics_locs, 1)
    for j = 1:size(mics_locs, 1)
        tau_test(i, j, :) = round((fs/c)*dot(repmat(mics_locs(i, :) - mics_locs(j, :), size(vMat, 1), 1) , vMat, 2));
    end
end

% first mic is reference
steering_test = zeros(size(mics_locs, 1), size(vMat, 1));
for i = 1:size(mics_locs, 1)
    steering_test(i, :) = round((fs/c)*dot(repmat(mics_locs(1, :) - mics_locs(i, :), size(vMat, 1), 1) , vMat, 2));
end

%%
for iter = 1:N
    test_point_norm = test_point(:, :, iter);
    
    hold on;plot3(test_point(:, 1, iter), test_point(:, 2, iter), test_point(:, 3, iter), 'rp', 'MarkerFaceColor', 'red')
    
    targetSteer = zeros(size(mics_locs, 1), size(test_point, 1));
    
    % test_point = vMat(15, :);
    % hold on;plot3(test_point(:, 1), test_point(:, 2), test_point(:, 3), 'rp', 'MarkerFaceColor', 'red')
    for j = 1:size(test_point, 1)
        for i = 1:size(mics_locs, 1)
            test_point_norm = test_point(j,:) / norm(test_point(j,:));
            targetSteer(i, j) = round((fs/c)*dot(mics_locs(1, :) - mics_locs(i, :) , test_point_norm, 2));
        end
    end
    signal1 = zeros(size(mics_locs,1),size(voice16K,1)+100);
    signal = signal1;
    for t = 1:length(targetIndx)
        for tt = 1:8
            signal1(tt,50+targetSteer(tt,t):size(voice16K,1)+targetSteer(tt,t)+49) = voice16K(:,t).';
        end
        signal = signal + P(t) * signal1;
    end
    
    FFT = fft(signal,[],2);
    FFTRep = repmat(FFT,8,1);
    FFTs = [];
    for ii = 1:8
        FFTs = [FFTs; FFTRep(ii:8:end,:)];
    end
    GCC = fftshift(ifft(FFTs .* conj(FFTRep),[],2),2);
    GCC(:,[1:end/2-99, end/2+100:end]) = [];
    
    [~,tauEst] = max(GCC.');
    tauEst = (tauEst - 100).';
    
    tau_test1 = reshape(tau_test,size(tau_test,1) * size(tau_test,2),[]);
    diff = 1./sum((tau_test1 - tauEst).^2).^0.5;
    [~, ind_max] = max(db(diff));
    
%     for j = 1:size(test_point, 1)
%         corr_test = repmat(targetSteer(:, j), 1, size(steering_test, 2)) - steering_test;
%         corr_test = sum(corr_test.^2).^0.5;
%         % corr_test = dot(repmat(test_steer, 1, size(steering_test, 2)) , steering_test, 1);
%         [~, ind_max(j)] = min(corr_test);
%     end
    
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







