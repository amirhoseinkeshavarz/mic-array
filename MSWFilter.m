%% MSW Filter
pqCrossCorrelationMswFilter = pqCrossCorrelation;
for iPQ = 1:micNum^2
    for iSample = delatTauPq(iPQ)+1:size(pqCrossCorrelation,2)-(delatTauPq(iPQ))
        pqCrossCorrelationMswFilter(iPQ,iSample) = max(pqCrossCorrelation(iPQ,iSample-delatTauPq(iPQ):iSample+delatTauPq(iPQ)));
    end
end

