%% directivity
gainMin = 0.1; % from table2 
alpha = 80; % from table2
beta = 100; % from table2
dMic = repmat(spacePointsCoarse(1,:),micNum,1);
u = spacePointsCoarse;

% equation 7 of paper
for i = 1:micNum
    for ii = 1:size(u,1)
        thetaUD(ii,i) = acosd(dot(u(ii,:),dMic(i,:))/norm(u(ii,:))/norm(dMic(i,:)));
    end
end

%equation 8 of paper
thetaUD = reshape(thetaUD,[],1);
gain = 1./ (1+exp(20/(beta - alpha) * (thetaUD - (alpha+beta)/2)));
gain = reshape(gain,size(spacePointsCoarse,1),micNum).';
gains = repmat(gain,micNum,1);
gain2 = [];
for i = 1:micNum
    gain2 = [gain2;gains(i:micNum:end,:)];
end

% equation 9 of paper
zeta = gains .* gain2;
zeta = zeta >= gainMin;



