function [x, mu, x1, P1, x2, P2] = one_step_imm2(x1, P1, x2, P2, z, F1, H1, Q1, F2, H2, Q2, R, mu, M)

    % 1 - Calculation of the mixing probabilities
    % cbar is the total probability, after interaction,
    % that the target is in state j. We use it as the
    % normalization constant.
    cbar = mu*M;

    omega = zeros(size(M));
    for j=1:length(M)
      for k=1:length(M)
        omega(j,k) = M(j, k)*mu(j)/cbar(k);
      end
    end

    % 2. Mixing

    % compute mixed initial conditions
    w1 = omega(:, 1);
    x1mix = zeros(length(x1), 1);
    x1mix = x1*w1(1) + x2*w1(2);
    P1mix = zeros(size(P1));
    y1 = x1-x1mix;

    w2 = omega(:, 2);
    x2mix = zeros(length(x1), 1);
    x2mix = x1*w2(1) + x2*w2(2);
    P2mix = zeros(size(P1));
    y2 = x2-x2mix;

    P1mix = w1(1)*(y1*y1' + P1) + w1(2)*(y2*y2' + P2);
    P2mix = w2(1)*(y1*y1' + P1) + w2(2)*(y2*y2' + P2);

    % 3 Mode-matched filtering

    % perform predict step using the mixed initial conditions
    [x1, P1] = kalman_predict(F1, Q1, 0, 0, x1mix, P1mix);

    % Update in kalman filter 1
    [x1, P1,residual1, S1, K] = kalman_update(H1, R, z, x1, P1);

    % perform predict step using the mixed initial conditions
    [x2, P2] = kalman_predict(F2, Q2, 0, 0, x2mix, P2mix);

    % Update in kalman filter 2
    [x2, P2, residual2, S2, K] = kalman_update(H2, R, z, x2, P2);

    % Likelihood
    L =  [mvnpdf(zeros(size(residual1)), residual1, S1), mvnpdf(zeros(size(residual2)), residual2, S2)];


     % 4 Mode probability update

     mu = cbar.*L;
     mu = mu/sum(mu, 2);

     % 5 Estimate and covariance combination.

     % compute mixed IMM state and covariance
     x = zeros(length(x1), 1);
     P = zeros(2,2);
     x = x1*mu(1) + x2*mu(2);
     y1 = x1-x;
     y2 = x2-x;
     P = mu(1)*(y1*y1' + P1) + mu(2)*(y2*y2' + P2);


end
