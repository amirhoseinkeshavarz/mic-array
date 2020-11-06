%% directivity
gainMin = 0.5;
alpha = 80;
beta = 100;
dMic = repmat(spacePoints(1,:),micNum,1);
u = spacePoints;
for i = 1:micNum
    for ii = 1:size(u,1)
        thetaUD(ii,i) = acosd(dot(u(ii,:),dMic(i,:))/norm(u(ii,:))/norm(dMic(i,:)));
    end
end

thetaUD = reshape(thetaUD,[],1);
gain = 1./ (1+exp(20/(beta - alpha) * (thetaUD - (alpha+beta)/2)));
gain = reshape(gain,size(spacePoints,1),micNum).';
gains = repmat(gain,micNum,1);
gain2 = [];
for i = 1:micNum
    gain2 = [gain2;gains(i:micNum:end,:)];
end

zeta = gains .* gain2;
zeta = zeta >= gainMin;



