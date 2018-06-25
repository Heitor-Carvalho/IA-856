function [x_est, likelyhood, residual_hist, trace_hist, Ppredict] = kalman_filter_wrapper(x0, P0, zv, F, H, Q, R)

  x = x0;                                     % Loading initial state mean
  P = P0;                                     % Loading initial cov matrix

  x_est = zeros(size(x,1), size(zv, 1));
  likelyhood = zeros(1, length(zv));
  trace_hist = zeros(1, length(zv));
  residual_hist = zeros(size(H,1), length(zv));

  for i=1:length(zv)
    z = zv(i, :)';

    % Predict step in kalman filter
    [x, P] = kalman_predict(F, Q, 0, 0, x, P);
    Ppredict = P;
    
    % Update
    [x, P,residual, S, K] = kalman_update(H, R, z, x, P);

    x_est(:, i) = x;

    % Likelihood
    likelyhood(i) =  mvnpdf(zeros(size(residual)), residual, S);

    % P trace
    trace_hist(i) = trace(P);

    % Residual
    residual_hist(:, i) = residual;
  end

end
