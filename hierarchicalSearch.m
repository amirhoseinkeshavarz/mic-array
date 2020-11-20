
for iCoarse = 1:size(TDMPsCoarse,2)
    if ~isempty(TDMPsCoarse(zeta(:,iCoarse),iCoarse))
        correlation(iCoarse) = 1./sum((TDMPsCoarse(zeta(:,iCoarse),iCoarse) - targetTDMP(zeta(:,iCoarse))).^2).^0.5;
    else
        correlation(iCoarse) = 0;
    end
end

%         correlation = 1./sum((TDMPs - targetTDMP).^2).^0.5;
[~ , indMaxCoarse(f,l)] = max(abs(correlation));

fineSearchNum = 10;
distancetoFinePoints = sum(abs(spacePointsFine - spacePointsCoarse(indMaxCoarse(f,l),:)).^2 ,2);
[distancetoFinePoints,indxSortDistance] = sort(distancetoFinePoints,'Ascend');
TDMPsFineTemp = TDMPsFine(:,indxSortDistance(1:fineSearchNum));

for iFine = 1:size(TDMPsFineTemp,2)
    if ~isempty(TDMPsCoarse(zeta(:,iFine),iFine))
        correlation(iFine) = 1./sum((TDMPsCoarse(zeta(:,iFine),iFine) - targetTDMP(zeta(:,iFine))).^2).^0.5;
    else
        correlation(iFine) = 0;
    end
end

%         correlation = 1./sum((TDMPs - targetTDMP).^2).^0.5;
[~ , indMaxFine(f,l)] = max(abs(correlation));
estimatedPosition = [spacePointsCoarse(indMaxFine(f,l),1),spacePointsCoarse(indMaxFine(f,l),2),spacePointsCoarse(indMaxFine(f,l),3)];

%     plot3(estimatedPosition(:, 1), estimatedPosition(:, 2), estimatedPosition(:, 3), 'bv', 'MarkerFaceColor', 'blue')
 