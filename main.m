clc
clear
close all

targetPositionIndx = [46 10];
targetNum = length(targetPositionIndx);
P = [1 1];
frameLength = 5000;
frameNumber = 16;

initialParams

Display

%% Time Difference of Microphone Pairs (TDMP)

% equation 12 of paper
TDMPsCoarse = zeros(size(micPosition, 1), size(micPosition, 1), size(spacePointsCoarse, 1));
for i = 1:size(micPosition, 1)
    for j = 1:size(micPosition, 1)
        TDMPsCoarse(i, j, :) = round((fs/c)*dot(repmat(micPosition(i, :) - micPosition(j, :)...
            , size(spacePointsCoarse, 1), 1) , spacePointsCoarse, 2));
    end
end
TDMPsCoarse = reshape(TDMPsCoarse , size(TDMPsCoarse,1) * size(TDMPsCoarse,2),[]);

TDMPsFine = zeros(size(micPosition, 1), size(micPosition, 1), size(spacePointsFine, 1));
for i = 1:size(micPosition, 1)
    for j = 1:size(micPosition, 1)
        TDMPsFine(i, j, :) = round((fs/c)*dot(repmat(micPosition(i, :) - micPosition(j, :)...
            , size(spacePointsFine, 1), 1) , spacePointsFine, 2));
    end
end
TDMPsFine = reshape(TDMPsFine , size(TDMPsFine,1) * size(TDMPsFine,2),[]);

%% Target Modelling

voiceImport

microphoneDirectivity

targetPosition = spacePointsCoarse(targetPositionIndx,:);

MSWCalib;

movement = 1*[0.09 -0.005 +0.01; -0.005 +0.09 +0.001];
% for t = 1:targetNum
%     targetMovement(t,:,:) = [linspace(0,movement(t,1),frameNumber).' linspace(0,movement(t,2),frameNumber).' linspace(0,movement(t,3),frameNumber).']';
% end
for f = 1:frameNumber
    % Received Signal
    hold on;
%     plot3(targetPosition(1, 1), targetPosition(1, 2), targetPosition(1, 3), 'mp', 'MarkerFaceColor', 'm') ; axis equal
%     plot3(targetPosition(2, 1), targetPosition(2, 2), targetPosition(2, 3), 'cp', 'MarkerFaceColor', 'c') ; axis equal
    targetTD = zeros(size(micPosition, 1), size(targetPosition, 1));
    for j = 1:size(targetPosition, 1)
        for i = 1:size(micPosition, 1)
            targetPositionNorm = targetPosition(j,:) / norm(targetPosition(j,:));
            targetTD(i, j) = round((fs/c)*dot(micPosition(1, :) - micPosition(i, :) , targetPositionNorm, 2));
        end
    end
    temp = zeros(size(micPosition,1),size(voice16K,1)+100);
    signal = temp;
    for t = 1:length(targetPositionIndx)
        for tt = 1:micNum
            temp(tt , 50 + targetTD(tt,t) : size(voice16K,1) + targetTD(tt,t) + 49) = gain(tt,targetPositionIndx(t)) * voice16K(:,f,t).';
        end
        signal = signal + P(t) * temp;
    end
    
    %% Target Localization 1
    
    SoundSourceLocalization
    
    targetPositionAll(f, :, :) = targetPosition;
    targetPosition = targetPosition + movement;
end
[indMaxFine Ed]

for ii = 1:targetNum
    estPosition(:,:,ii) = spacePointsFine(indMaxFine(:,ii),:);
end
nTarget = size(estPosition, 3);
nMeasurement = size(estPosition, 3);

%% test target modeling
% for iFrame = 1:size(targetPositionAll, 1)
%     for iTarget = 1:size(targetPositionAll, 2)
%     estPosition(iFrame, :, iTarget) = targetPositionAll(iFrame, iTarget, :) / ...
%                                 norm(squeeze(targetPositionAll(iFrame, iTarget, :)));
%     end
% end
Ed = [0.2*ones(size(targetPositionAll, 1), 1), 0.2*ones(size(targetPositionAll, 1), 1)];
% Ed(5, 1) = 0.01;
% estPosition(5, :, 1) = 2*estPosition(5, :, 1);


%% tracking
% parameters
delta_T = 1;
F = eye(6);
for i = 4:6
    F(i-3, i) = delta_T;
end

sigma2_Q = 0.000009;
Q = zeros(6);
for i = 4:6
    Q(i, i) =  sigma2_Q;
end

H = zeros(3, 6);
for i = 1:3
    H(i, i) = 1;
end

sigma2_R = 0.0015;
R = sigma2_R*eye(3);

sigma2_R_Prob = 0.003;
R_prob = sigma2_R_Prob*eye(3);

% initial state
state_initial = [squeeze(estPosition(1, :, :)).', movement].';
for iTarget = 1:nTarget
    P_initial(:, :, iTarget) = 0.1*eye(6); % what should be its initial value?
end

state_posterior = state_initial;
P_posterior = P_initial;
nTarget; % this should be set during tracking
nMeasurement; % this should be set during tracking
for iter = 1:frameNumber
%     plot3(squeeze(estPosition(iter, 1, 1)), squeeze(estPosition(iter, 2, 1)), squeeze(estPosition(iter, 3, 1)), 'mv', 'MarkerFaceColor', 'm')
%     plot3(squeeze(estPosition(iter, 1, 2)), squeeze(estPosition(iter, 2, 2)), squeeze(estPosition(iter, 3, 2)), 'cv', 'MarkerFaceColor', 'c')

%     plot3(state_posterior(1, 1), state_posterior(2, 1), state_posterior(3, 1), 'm-o')
%     plot3(state_posterior(1, 2), state_posterior(2, 2), state_posterior(3, 2), 'c-o')
    state_posteriorAll(:, :, iter) = state_posterior;

    % Prediction(Step A)
    state_predicted = F*state_posterior;
    for iTarget = 1:nTarget
        P_predicted(:, :, iTarget) = F*P_posterior(:, :, iTarget)*F' + Q;
    end
    
    % Normalization(Step B)
    d = state_predicted(1:3, :);
    s = state_predicted(4:6, :);
    
    for j = 1:nTarget
        d_norm(:, j) = d(:, j) / norm(d(:, j));
        s_norm(:, j) = s(:, j) - d(:, j) * (s(:, j).' * d(:, j) / norm(d(:, j))^2);
    end
    state_predicted = [d_norm; s_norm];
    
    % Assignment(Step C)
    assignment_alphabet = [-2, -1, 1:nTarget];
    assignment_alphabet = [assignment_alphabet, assignment_alphabet];
    Fg = perms(assignment_alphabet);
    Fg = Fg(:, 1:2);
    Fg = unique(Fg, 'rows');
    
    % Likelihood(Step D)
    
    mio_i = H*state_predicted;
    for iTarget = 1:nTarget
        Sigma_i(:, :, iTarget) = H*P_predicted(:, :, iTarget)*H';
    end
    
    mio_v = squeeze(estPosition(iter, :, :));
    for iTarget = 1:nTarget
        Sigma_v(:, :, iTarget) = R;
    end
    
    for iTrack = 1:size(mio_i, 2)
        for vMeasurement = 1:size(mio_v, 2)
            Sigma_iv(:, :, iTrack, vMeasurement) = (Sigma_i(:, :, iTrack)^-1 + Sigma_v(:, :, vMeasurement)^-1)^-1;
            mio_iv(:, iTrack, vMeasurement) = Sigma_iv(:, :, iTrack, vMeasurement)*...
                (Sigma_i(:, :, iTrack)^-1*mio_i(:, iTrack) + Sigma_v(:, :, vMeasurement)^-1*mio_v(:, vMeasurement));
            
            C1_iv(iTrack, vMeasurement) = log(det(Sigma_iv(:, :, iTrack, vMeasurement)))...
                                        - log(8*(pi^3)*det(Sigma_i(:, :, iTrack))*det(Sigma_v(:, :, vMeasurement)));
            C2_iv(iTrack, vMeasurement) = mio_iv(:, iTrack, vMeasurement)' * Sigma_iv(:, :, iTrack, vMeasurement)^-1 *...
                                          mio_iv(:, iTrack, vMeasurement);
            C3_iv(iTrack, vMeasurement) = mio_i(:, iTrack)' * Sigma_i(:, :, iTrack)^-1 * mio_i(:, iTrack);
            C4_iv(iTrack, vMeasurement) = mio_v(:, vMeasurement)' * Sigma_v(:, :, vMeasurement)^-1 * mio_v(:, vMeasurement);
            
            omega_iv(iTrack, vMeasurement) = exp(0.5*(C1_iv(iTrack, vMeasurement) + C2_iv(iTrack, vMeasurement) - ...
                                                 C3_iv(iTrack, vMeasurement) - C4_iv(iTrack, vMeasurement)));
                                             
            Prob_lambdav_Ci(iTrack, vMeasurement) = omega_iv(iTrack, vMeasurement); 
        end
    end
 
    mu_Inactive = 0.1; 
    sigma_Inactive = 0.0025;
%     Prob_LAMBDA_Inactive = @(LAMBDA) (1/sqrt(2*pi*sigma_Inactive^2))*(exp(-(LAMBDA - mu_Inactive)^2/ (2*sigma_Inactive^2)));
    Prob_LAMBDA_Inactive = @(LAMBDA) pdf(makedist('Normal','mu',mu_Inactive,'sigma',sigma_Inactive), LAMBDA);

    mu_Active = 0.2; 
    sigma_Active = 0.0025;
    Prob_LAMBDA_Active = @(LAMBDA)  pdf('Normal', LAMBDA, mu_Active, sigma_Active);
    
    Prob_LAMBDA_Diffused = 0.005;  % reconsider
    
    Prob_ksiv_fgv = zeros(size(Fg, 1), nMeasurement);
    for g = 1:size(Fg, 1)
        for v = 1:nMeasurement
            if Fg(g, v) == -2
                Prob_ksiv_fgv(g, v) = Prob_LAMBDA_Inactive(Ed(iter, v)) * Prob_LAMBDA_Diffused;
            end
            if Fg(g, v) == -1
                Prob_ksiv_fgv(g, v) = Prob_LAMBDA_Active(Ed(iter, v)) * Prob_LAMBDA_Diffused;
            end
            if Fg(g, v) >= 1
                i = Fg(g, v);
                Prob_ksiv_fgv(g, v) = Prob_LAMBDA_Active(Ed(iter, v)) * Prob_lambdav_Ci(i, v);
            end
        end
    end
    
    Prob_KSI_Fg = prod(Prob_ksiv_fgv, 2) ;
    
    % Prior probability(Step E)
    P_new = 0.1; % reconsider
    P_false = 0.1; % reconsider
    P_track = 0.8; % reconsider
    
    for g = 1:size(Fg, 1)
        for v = 1:nMeasurement
            if Fg(g, v) == -2
                Prob_fg(g, v) = P_false;
            end
            if Fg(g, v) == -1
                Prob_fg(g, v) = P_new;
            end
            if Fg(g, v) >= 1
                Prob_fg(g, v) = P_track;
            end
        end
    end
    
    Prob_Fg = prod(Prob_fg, 2);
    
    % Poseterior Prabability(Step F)
    
    Bayes_denominator = sum(Prob_KSI_Fg.*Prob_Fg);
    Prob_Fg_KSI = Prob_KSI_Fg .* Prob_Fg / Bayes_denominator;
    
    for vMeasurement = 1:nMeasurement
        for iTrack = 1:nTarget
            ind_g = find(Fg(:, vMeasurement) == iTrack);
            Prob_i_KSIv(iTrack, vMeasurement) = sum(Prob_Fg_KSI(ind_g));
        end
        ind_new_g = find(Fg(:, vMeasurement) == -1);
        Prob_new_KSIv(vMeasurement) = sum(Prob_Fg_KSI(ind_new_g));
    end
    
    for iTrack = 1:nTarget
        ind_g = logical(sum(Fg == iTrack, 2) >= 1);
        Prob_i_KSI(iTrack) = sum(Prob_Fg_KSI(ind_g));
    end
    
    
    % Adding and removing sources (Step G)
    
    
    % Update (Step H)
    for iTrack = 1:nTarget
        K(:, :, iTrack) = P_predicted(:, :, iTrack) * H' * ( H*P_predicted(:, :, iTrack)*H' + R )^-1; % Kalman gain
        [~, vhat(iTrack)] = max(Prob_i_KSIv(iTrack, :));
        
        
        state_posterior(:, iTrack) = state_predicted(:, iTrack) + Prob_i_KSI(iTrack) * K(:, :, iTrack) * ( estPosition(iter, :, vhat(iTrack)).' - H*state_predicted(:, iTrack) );
        P_posterior(:, :, iTrack) = P_predicted(:, :, iTrack) - Prob_i_KSI(iTrack) * K(:, :, iTrack) * H * P_predicted(:, :, iTrack);
    end

    %     plot3(estimated_locs(:, 1), estimated_locs(:, 2), estimated_locs(:, 3), 'bv', 'MarkerFaceColor', 'blue')
    
end

plot3(squeeze(state_posteriorAll(1, 1, :)), squeeze(state_posteriorAll(2, 1, :)), ...
    squeeze(state_posteriorAll(3, 1, :)), 'm-o', 'MarkerFaceColor', [0.25 0.34 0.28], ...
    'MarkerSize', 8, 'DisplayName', 'Track1')
                                        
plot3(squeeze(state_posteriorAll(1, 2, :)), squeeze(state_posteriorAll(2, 2, :)), ...
    squeeze(state_posteriorAll(3, 2, :)), 'c-o', 'MarkerFaceColor', [0.85 0.3 0.1], ...
    'MarkerSize', 8, 'DisplayName', 'Track2')


plot3(squeeze(estPosition(:, 1, 1)), squeeze(estPosition(:, 2, 1)), squeeze(estPosition(:, 3, 1)), 'mv', 'MarkerFaceColor', 'm', 'DisplayName', 'estimated position 1')
plot3(squeeze(estPosition(:, 1, 2)), squeeze(estPosition(:, 2, 2)), squeeze(estPosition(:, 3, 2)), 'cv', 'MarkerFaceColor', 'c', 'DisplayName', 'estimated position 2')

plot3(squeeze(targetPositionAll(:, 1, 1)), squeeze(targetPositionAll(:, 1, 2)), squeeze(targetPositionAll(:, 1, 3)), 'm-p', 'MarkerFaceColor', 'm', 'DisplayName', 'target location 1') ; axis equal
plot3(squeeze(targetPositionAll(:, 2, 1)), squeeze(targetPositionAll(:, 2, 2)), squeeze(targetPositionAll(:, 2, 3)), 'c-p', 'MarkerFaceColor', 'c', 'DisplayName', 'target location 2') ;

legend('show')



