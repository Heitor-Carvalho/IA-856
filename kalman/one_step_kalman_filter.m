function [x, P, likelyhood] = one_step_kalman_filter(x, P, z, F, H, Q, R)

  % Predict step in kalman filter
  [x, P] = kalman_predict(F, Q, 0, 0, x, P);

  % Update
  [x, P,residual, S, K] = kalman_update(H, R, z, x, P);

  % Likelihood
  likelyhood =  mvnpdf(zeros(size(residual)), residual, S);

end
