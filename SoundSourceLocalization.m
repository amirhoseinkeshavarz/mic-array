for l = 1:targetNum
    % equation 3 from paper
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
    pqCrossCorrelation = fftshift(ifft(FFTs .* conj(FFTRep)./(FFTNorm),[],2),2);
    
    % equation 22 from paper
    MSWFilter;
    
    pqCrossCorrelationMswFilter(:,[1:end/2-99, end/2+100:end]) = [];
    [~,targetTDMP] = max(pqCrossCorrelationMswFilter.');
    targetTDMP = (targetTDMP - 100).';
    
    hierarchicalSearch;
    
    for ll = 1:micNum
        signal1(ll,:) = circshift(signal(ll,:),TDMPsFine(ll,indMaxFine(f,l)));
    end
    signal11 = mean(signal1,1);
    for ll = 1:micNum
        signal22(ll,:) = circshift(signal11,-TDMPsFine(ll,indMaxFine(f,l)));
    end
    signal = signal - signal22;
    Ed(f,l) = sqrt(sum(abs(signal11.^2)));
end

clear signal
