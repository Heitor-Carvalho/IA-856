function [x, P] = kalmna_pedict(F, Q, B, u, x, P)
  % Predict step in kalman filter
  if(B == 0)
    x = F*x;            % Propagation of the mean though state equation
  else
    x = F*x + B*u
  end

  P = F*P*F' + Q;       % Propagation of covariance matrix though state equation
end
