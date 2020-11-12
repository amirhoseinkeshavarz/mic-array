clc
clear
close all

targetPositionIndx = [46 10];
targetNum = length(targetPositionIndx);
P = [1 1];
frameLength = 256;
frameNumber = 300;

initialParams

Display

%% Time Difference of Microphone Pairs (TDMP)

TDMPs = zeros(size(micPosition, 1), size(micPosition, 1), size(spacePoints, 1));
for i = 1:size(micPosition, 1)
    for j = 1:size(micPosition, 1)
        TDMPs(i, j, :) = round((fs/c)*dot(repmat(micPosition(i, :) - micPosition(j, :), size(spacePoints, 1), 1) , spacePoints, 2));
    end
end
TDMPs = reshape(TDMPs , size(TDMPs,1) * size(TDMPs,2),[]);

%% Target Modelling

voiceImport
microphoneDirectivity
targetPosition = spacePoints(targetPositionIndx,:);
movement = [0.001 +0.002 +0.01; -0.01 +0.02 +0.01];
for t = 1:targetNum
    targetMovement(:,:,t) = [linspace(0,movement(t,1),frameNumber).' linspace(0,movement(t,2),frameNumber).' linspace(0,movement(t,3),frameNumber).'];
end


    MSWCalib;


for f = 1:frameNumber
    % Received Signal
    hold on;plot3(targetPosition(:, 1), targetPosition(:, 2), targetPosition(:, 3), 'rp', 'MarkerFaceColor', 'red')
    targetTD = zeros(size(micPosition, 1), size(targetPosition, 1));
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
            temp(tt , 50 + targetTD(tt,t) : size(voice16K,1) + targetTD(tt,t) + 49) = gain(tt,targetPositionIndx(t)) * voice16K(:,f,t).';
        end
        signal = signal + P(t) * temp;
    end
    
    
    %% Target Localization 1
    
    SoundSourceLocalization
    targetPosition = targetPosition + squeeze(targetMovement(f,:,:)).';
end
[indMax Ed]

for ii = 1:targetNum
    estPosition(:,:,ii) = spacePoints(indMax(:,ii),:);
end

%}
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



