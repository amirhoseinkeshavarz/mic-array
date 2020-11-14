micNum = 8;
fs = 16000;
c = 343 + 5 * randn(1) % m/s speed of sound
L = 3;
spacePointNum = 10*4^L + 2;
%% microphones location
theta = (0 : 2*pi/micNum : 2*pi - 2*pi/micNum)';
arrayRadius = 0.5/2;
micPosition = arrayRadius * [cos(theta), sin(theta), zeros(length(theta), 1)];
% micPosition = 
%% creating test points
[spacePoints, fMat] = spheretri(spacePointNum);
zMinus = spacePoints(:,3)<0;
spacePoints(zMinus,:) = []; 