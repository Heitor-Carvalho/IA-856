function [x_est, mu_est] = imm3(x0, P0, zv, F1, H1, Q1, F2, H2, Q2, F3, H3, Q3, R, mu, M)

  x1 = x0;                                     % Loading initial state mean 1
  P1 = P0;                                     % Loading initial cov matrix 1
  x2 = x0;                                     % Loading initial state mean 2
  P2 = P0;                                     % Loading initial cov matrix 2
  x3 = x0;                                     % Loading initial state mean 2
  P3 = P0;                                     % Loading initial cov matrix 2

  x_est = zeros(size(x0,1), size(zv, 1));
  mu_est = zeros(size(mu,2), size(zv, 1));

  for i=1:length(zv)
    z = zv(i, :)';

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
    x1mix = zeros(length(x0), 1);
    x1mix = x1*w1(1) + x2*w1(2) + x3*w1(3);
    P1mix = zeros(size(P0));
    y1 = x1-x1mix;

    w2 = omega(:, 2);
    x2mix = zeros(length(x0), 1);
    x2mix = x1*w2(1) + x2*w2(2) + x3*w2(3);
    P2mix = zeros(size(P0));
    y2 = x2-x2mix;

    w3 = omega(:, 3);
    x3mix = zeros(length(x0), 1);
    x3mix = x1*w3(1) + x2*w3(2) + x3*w3(3);
    P3mix = zeros(size(P0));
    y3 = x3-x3mix;

    P1mix = w1(1)*(y1*y1' + P1) + w1(2)*(y2*y2' + P2)+ w1(3)*(y3*y3' + P3);
    P2mix = w2(1)*(y1*y1' + P1) + w2(2)*(y2*y2' + P2)+ w2(3)*(y3*y3' + P3);
    P3mix = w3(1)*(y1*y1' + P1) + w3(2)*(y2*y2' + P2)+ w3(3)*(y3*y3' + P3);

    % 3 Mode-matched filtering

    % perform predict step using the mixed initial conditions
    [x1, P1] = kalman_predict(F1, Q1, 0, 0, x1mix, P1mix);

    % Update in kalman filter 1
    [x1, P1,residual1, S1, K] = kalman_update(H1, R, z, x1, P1);

    % perform predict step using the mixed initial conditions
    [x2, P2] = kalman_predict(F2, Q2, 0, 0, x2mix, P2mix);

    % Update in kalman filter 2
    [x2, P2, residual2, S2, K] = kalman_update(H2, R, z, x2, P2);

    % perform predict step using the mixed initial conditions
    [x3, P3] = kalman_predict(F3, Q3, 0, 0, x3mix, P3mix);

    % Update in kalman filter 2
    [x3, P3, residual3, S3, K] = kalman_update(H3, R, z, x3, P3);

    % Likelihood
    L =  [mvnpdf(zeros(size(residual1)), residual1, S1), mvnpdf(zeros(size(residual2)), residual2, S2), mvnpdf(zeros(size(residual3)), residual3, S3)];


     % 4 Mode probability update

     mu = cbar.*L;
     mu = mu/sum(mu, 2);

     % 5 Estimate and covariance combination.

     % compute mixed IMM state and covariance
     x = zeros(length(x0), 1);
     P = zeros(2,2);
     x = x1*mu(1) + x2*mu(2);
     y1 = x1-x;
     y2 = x2-x;
     y3 = x3-x;
     P = mu(1)*(y1*y1' + P1) + mu(2)*(y2*y2' + P2) + mu(2)*(y3*y3' + P3);

    x_est(:, i) = x;
    mu_est(:, i) = mu;
  end

end

