clc
clear
close all

targetPositionIndx = [250];
targetNum = length(targetPositionIndx);
P = [1 1 1];
iterNum = 1;
initialParams

voiceImport

%% Target Modelling
% test_vel = [0.5*(rand(test_points_NUM, 2) - 0.5) rand(test_points_NUM, 1)];

targetPosition = spacePoints(targetPositionIndx,:);

Display

%% Time Difference of Microphone Pairs (TDMP)

TDMPs = zeros(size(micPosition, 1), size(micPosition, 1), size(spacePoints, 1));
for i = 1:size(micPosition, 1)
    for j = 1:size(micPosition, 1)
        TDMPs(i, j, :) = round((fs/c)*dot(repmat(micPosition(i, :) - micPosition(j, :), size(spacePoints, 1), 1) , spacePoints, 2));
    end
end
TDMPs = reshape(TDMPs , size(TDMPs,1) * size(TDMPs,2),[]);

%%
for iter = 1:iterNum
    targetPositionNorm = targetPosition(:, :, iter);
    
    hold on;plot3(targetPosition(:, 1, iter), targetPosition(:, 2, iter), targetPosition(:, 3, iter), 'rp', 'MarkerFaceColor', 'red')
    
    targetTD = zeros(size(micPosition, 1), size(targetPosition, 1));
    
    % test_point = vMat(15, :);
    % hold on;plot3(test_point(:, 1), test_point(:, 2), test_point(:, 3), 'rp', 'MarkerFaceColor', 'red')
    for j = 1:size(targetPosition, 1)
        for i = 1:size(micPosition, 1)
            targetPositionNorm = targetPosition(j,:) / norm(targetPosition(j,:));
            targetTD(i, j) = round((fs/c)*dot(micPosition(1, :) - micPosition(i, :) , targetPositionNorm, 2));
        end
    end
    temp = zeros(size(micPosition,1),size(voice16K,1)+100);
    signal = temp;
    for t = 1:length(targetPositionIndx)
        for tt = 1:micNum
            temp(tt , 50 + targetTD(tt,t) : size(voice16K,1) + targetTD(tt,t) + 49) = voice16K(:,t).';
        end
        signal = signal + P(t) * temp;
    end
    
    %% Target Localization
    for l = 1:targetNum
        FFT = fft(signal,[],2);
        FFTRep = repmat(FFT,micNum,1);
        FFTRepNorm = sqrt(sum(FFTRep.^2,2));
        FFTs = [];
        FFTsNorm = [];
        for ii = 1:micNum
            FFTs = [FFTs; FFTRep(ii:micNum:end,:)];
            FFTsNorm = [FFTsNorm ; FFTRepNorm(ii:micNum:end)];
        end
        FFTNorm = FFTsNorm .* FFTRepNorm;
        GCC = fftshift(ifft(FFTs .* conj(FFTRep)./(FFTNorm),[],2),2);
        GCC(:,[1:end/2-99, end/2+100:end]) = [];
        
        [~,targetTDMP] = max(GCC.');
        targetTDMP = (targetTDMP - 100).';
        
        correlation = 1./sum((TDMPs - targetTDMP).^2).^0.5;
        [~, indMax] = max(db(correlation));
        estimatedPosition = [x(indMax),y(indMax),z(indMax)];
        plot3(estimatedPosition(:, 1), estimatedPosition(:, 2), estimatedPosition(:, 3), 'bv', 'MarkerFaceColor', 'blue')
        for ll = 1:8
            signal1(ll,:) = circshift(signal(ll,:),TDMPs(ll,indMax));
        end
        signal11 = mean(signal1,1);
        for ll = 1:8
            signal22(ll,:) = circshift(signal11,-TDMPs(ll,indMax));
        end
        signal = signal - signal22;
    end
end

%% tracking
%{
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
%}






