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
    [~ , indMax(f,l)] = max(abs(correlation));
    estimatedPosition = [x(indMax(f,l)),y(indMax(f,l)),z(indMax(f,l))];
%     plot3(estimatedPosition(:, 1), estimatedPosition(:, 2), estimatedPosition(:, 3), 'bv', 'MarkerFaceColor', 'blue')
    
    for ll = 1:micNum
        signal1(ll,:) = circshift(signal(ll,:),TDMPs(ll,indMax(f,l)));
    end
    signal11 = mean(signal1,1);
    for ll = 1:micNum
        signal22(ll,:) = circshift(signal11,-TDMPs(ll,indMax(f,l)));
    end
    signal = signal - signal22;
    Ed(f,l) = sqrt(sum(abs(signal11.^2)));
end

clear signal
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