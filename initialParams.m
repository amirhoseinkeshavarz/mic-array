micNum = 8;
fs = 16000;
L = 3;
cMin = 0.3;  % table2 of paper
micVariance = 1e-6; % table2 of paper
muC = 343; % table2 of paper 
sigmaC = 5; % table2 of paper
c = muC + sigmaC * randn(1); % m/s speed of sound

%% microphones location
theta = (0 : 2*pi/micNum : 2*pi - 2*pi/micNum)';
arrayRadius = 0.5/2;
micPosition = arrayRadius * [cos(theta), sin(theta), zeros(length(theta), 1)];
micPosition = micPosition + sqrt(micVariance) * randn(size(micPosition));

%% creating test points
spacePointCoarseNum = 10*4^L + 2;
[spacePointsCoarse, fMatCoarse] = spheretri(spacePointCoarseNum);
zMinusCoarse = spacePointsCoarse(:,3)<0;
spacePointsCoarse(zMinusCoarse,:) = []; 

spacePointFineNum = 10*4^(L+1) + 2;
[spacePointsFine, fMatFine] = spheretri(spacePointFineNum);
zMinusFine = spacePointsFine(:,3)<0;
spacePointsFine(zMinusFine,:) = []; 
