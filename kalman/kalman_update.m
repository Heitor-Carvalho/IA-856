function [x, P,residual, S, K] = kalman_update(H, R, z, x, P)
  % Update step in kalman filter
  residual = z - H*x;                   % Measurment residual

  PHT = P*H';                           % Common term used

  S = H*PHT + R;                        % Measurment covariance

  K = PHT*inv(S);                       % Filter gain

  x = x + K*residual;                   % State mean update

  I = eye(length(x));
  I_KH = I - K*H;
  P = I_KH*P*I_KH' + K*R*K';           % Covariance update

end
