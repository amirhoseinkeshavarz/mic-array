function [tauProbabilitySum,vuProbability] = vuProbabilityFunc(u,micPosition,pq,sigmaC,muC,micNum,fs,delatTauPq)
for iU = 1:size(u,1)
    octagon;
    for iOctagon = 1:size(octagonXyz,1)
        cntPq = 1;
        for iP = 1:micNum
            for iQ = 1:micNum
                
                pq.deltaTau = delatTauPq(iP * iQ); % initial value is zero (step2 of Algorithm1)
                
                % equation 21 of paper
                pq.tauEstimated = round(fs/muC * dot((micPosition(iP,:) - micPosition(iQ,:)) , octagonXyz(iOctagon,:)));
                
                % equation 23 of paper
                alpha = pq.tauEstimated - pq.deltaTau - 0.5;
                
                % equation 24 of paper
                beta = pq.tauEstimated + pq.deltaTau + 0.5;
                
                % equation 19 of paper
                pq.tauMean = fs/muC * dot((micPosition(iP,:) - micPosition(iQ,:)) ,octagonXyz(iOctagon,:));
                
                % equation 20 of paper
                pq.standardDeviation = fs/muC * sqrt(octagonXyz(iOctagon,:) * (pq.crossStd * 2) * octagonXyz(iOctagon,:)' +...
                    dot((micPosition(iP,:) - micPosition(iQ,:)) , octagonXyz(iOctagon,:))^2 * sigmaC^2/muC^2 );
                
                % equation 25 of paper
                fun = @(tau) 1/sqrt(2*pi*pq.standardDeviation^2)* exp(-(tau-pq.tauMean).^2/2/pq.standardDeviation^2);
                tauProbability(cntPq,iOctagon,iU) = integral(fun,alpha,beta);
                
                cntPq = cntPq + 1;
            end
        end
    end
    iU
end
% equation 27 of paper
tauProbabilitySum = sum(sum(tauProbability,3),2)/nNeighbour/size(u,1); 

tauP = [];

% equation 28 of paper
for iP = 1 : micNum
    tauP = [tauP; tauProbability((iP-1) * micNum + iP +1:iP * micNum,:,:)];
end
tauP = tauProbability;
tauP(1:8:end,:,:) = [];
vuProbability = squeeze(sum(tauP,1))/micNum/(micNum-1); 


end