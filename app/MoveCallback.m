function MoveCallback(object, eventdata)
  global F
  global H
  global Q
  global R
  global x
  global P
  global pos_noise_hist
  global pos_hist
  global x_hist
  global lk_hist
  global sigma_r

  C =  get(0,'PointerLocation');
  z = C';

  pos_hist = circshift(pos_hist, 1);
  pos_hist(1, :) = z;

  z_noise = z + sigma_r*randn(size(z));
  pos_noise_hist = circshift(pos_noise_hist, 1);
  pos_noise_hist(1, :) = z_noise;

  [x, P, lk] = one_step_kalman_filter(x, P, z_noise, F, H, Q, R);

  lk_hist = circshift(lk_hist, 1);
  lk_hist(1, :) = lk;


  x_hist = circshift(x_hist, 1);
  x_hist(1, :) = x;

  figure(2)
  clf
  plot(lk_hist, 'o--')
  ylim([0.1 10]*1e-4)
  grid on

  figure(1)
  clf
  plot(pos_noise_hist(:, 1), pos_noise_hist(:, 2), 'ob')
  hold on
  plot(pos_hist(:, 1), pos_hist(:, 2), 'k')
  plot(x_hist(:, 1), x_hist(:, 3), '.--r')
  axis([0 1 0 1]*1000/2)
  grid on


end
