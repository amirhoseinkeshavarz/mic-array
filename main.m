clc
clear
close all

targetPositionIndx = [10 46];
targetNum = length(targetPositionIndx);
P = [1 1];
frameLength = 250000;

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
        temp(tt , 50 + targetTD(tt,t) : size(voice16K,1) + targetTD(tt,t) + 49) = gain(tt,targetPositionIndx(t)) * voice16K(:,t).';
    end
    signal = signal + P(t) * temp;
end

%% Target Localization 2
%{
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
    FFTCorr = FFTs .* conj(FFTRep)./(FFTNorm);
    for nn = 1:size(TDMPs,2)
        Rij1 = FFTCorr .* exp(1i* 2*pi * TDMPs(:,nn)*(0:size(FFTs,2)-1)/size(FFTs,2));
        Rij(:,nn) = sum(Rij1,2);
    end
    
    Ed = sum(abs(Rij));
    [targetEnergy(l),indMax(l)] = max(Ed);
    
    TDMPs(:,IndMax(l)) = zeros(size(TDMPs,1),1);
end
indMax
%}
%% Target Localization 1
%
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
    for kk = 1:size(TDMPs,2)
        if ~isempty(TDMPs(zeta(:,kk),kk))
            correlation(kk) = 1./sum((TDMPs(zeta(:,kk),kk) - targetTDMP(zeta(:,kk))).^2).^0.5;
        else
            correlation(kk) = 0;
        end
    end

%         correlation = 1./sum((TDMPs - targetTDMP).^2).^0.5;
    [~ , indMax(l)] = max(abs(correlation));
    estimatedPosition = [x(indMax(l)),y(indMax(l)),z(indMax(l))];
    plot3(estimatedPosition(:, 1), estimatedPosition(:, 2), estimatedPosition(:, 3), 'bv', 'MarkerFaceColor', 'blue')
    
    for ll = 1:8
        signal1(ll,:) = circshift(signal(ll,:),TDMPs(ll,indMax(l)));
    end
    signal11 = mean(signal1,1);
    for ll = 1:8
        signal22(ll,:) = circshift(signal11,-TDMPs(ll,indMax(l)));
    end
    signal = signal - signal22;
    Ed(l) = sqrt(sum(abs(signal11.^2)));
end
indMax
Ed
%}

