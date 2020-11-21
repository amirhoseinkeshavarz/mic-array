nNeighbour = 8;
% xyzDistance = 0.1 * randn(3,nNeighbour);
% octagonXyz = u(iU,:) + xyzDistance.';

radius = 0.1;
theta = 2*pi/nNeighbour : 2*pi/nNeighbour :2*pi;
v=null(u(iU,:));
octagonXyz = (repmat(u(iU,:)',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta)))';

for iOcta = 1:size(octagonXyz,1)
    octagonXyz(iOcta,:) = octagonXyz(iOcta,:)/norm(octagonXyz(iOcta,:));
end
