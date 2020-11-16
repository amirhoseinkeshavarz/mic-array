cMin = 0.3;  % table2 of paper
micVariance = 1e-6; % table2 of paper
muC = 343; % table2 of paper 
sigmaC = 5; % table2 of paper
pq.crossStd = micVariance * eye(3);


delatTauPq = zeros(micNum^2,1);
vuProbability = zeros(micNum,size(u,1));
pqStar = [];
while min(min(vuProbability)) < cMin
delatTauPq(pqStar) = delatTauPq(pqStar) + 1;

for iU = 1:size(u,1)
    octagon;
    for iOctagon = 1:size(octagonXyz,1)
        cntPq = 1;
        for iP = 1:micNum
            for iQ = 1:micNum
                
                pq.deltaTau = delatTauPq(iP * iQ);
                pq.tauEstimated = fs/muC * dot((micPosition(iP,:) - micPosition(iQ,:)) , u(iU,:));
                
                alpha = pq.tauEstimated - pq.deltaTau - 0.5;
                beta = pq.tauEstimated + pq.deltaTau + 0.5;
                
                pq.tauMean = fs/muC * dot((micPosition(iP,:) - micPosition(iQ,:)) , u(iU,:));
                pq.standardDeviation = fs/muC * sqrt(u(iU,:) * (pq.crossStd * 2) * u(iU,:)' +...
                    dot((micPosition(iP,:) - micPosition(iQ,:)) , u(iU,:))^2 * sigmaC^2/muC^2 );
                
                fun = @(tau) 1/sqrt(2*pi*pq.standardDeviation^2)* exp(-(tau-pq.tauMean).^2/2/pq.standardDeviation^2);
                tauProbability(cntPq,iOctagon,iU) = integral(fun,alpha,beta);
                cntPq = cntPq + 1;
            end
        end
    end
    iU
end

tauProbabilitySum = sum(sum(tauProbability,3),2)/nNeighbour/size(u,1); %% equation 27 of paper

tauP = [];
for iP = 1 : micNum
    tauP = [tauP; tauProbability((iP-1) * micNum + iP +1:iP * micNum,:,:)];
end
vuProbability = squeeze(sum(tauP,1))/micNum/(micNum-1); %% equation 28 of paper
%%
[pqValue,pqStar] = min(tauProbabilitySum);

end


