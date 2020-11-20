
pq.crossStd = micVariance * eye(3);


delatTauPq = zeros(micNum^2,1);
vuProbability = zeros(micNum,size(u,1));
pqStar = [];
[tauProbabilitySum,vuProbability] = vuProbabilityFunc(u,micPosition,pq,sigmaC,muC,micNum,fs,delatTauPq);
while min(min(vuProbability)) < cMin
delatTauPq(pqStar) = delatTauPq(pqStar) + 1;
[pqValue,pqStar] = min(tauProbabilitySum);
[tauProbabilitySum,vuProbability] = vuProbabilityFunc(u,micPosition,pq,sigmaC,muC,micNum,fs,delatTauPq);

end


