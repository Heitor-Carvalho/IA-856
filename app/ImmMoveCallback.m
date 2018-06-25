function ImmMoveCallback(object, eventdata)
  global F1
  global H1
  global Q1
  global x1
  global P1
  global F2
  global H2
  global Q2
  global x2
  global P2
  global pos_noise_hist
  global pos_hist
  global x_hist
  global mu_hist
  global R
  global sigma_r
  global mu
  global M

  C =  get(0,'PointerLocation');
  z = C';
  
  pos_hist = circshift(pos_hist, 1);
  pos_hist(1, :) = z;

  z_noise = z + sigma_r*randn(size(z));
  pos_noise_hist = circshift(pos_noise_hist, 1);
  pos_noise_hist(1, :) = z_noise;

  [x, mu, x1, P1, x2, P2] = one_step_imm2(x1, P1, x2, P2, z_noise, F1, H1, Q1, F2, H2, Q2, R, mu, M);

  x_hist = circshift(x_hist, 1);
  x_hist(1, :) = x;

  mu_hist = circshift(mu_hist, 1);
  mu_hist(1, :) = mu;

  figure(2)
  clf
  plot(mu_hist(:, 1), 'ro--')
  hold on
  plot(mu_hist(:, 2), 'bo--')
  ylim([0 1])
  grid on

  figure(1)
  clf
  plot(pos_noise_hist(:, 1), pos_noise_hist(:, 2), 'ob')
  hold on
  plot(pos_hist(:, 1), pos_hist(:, 2), 'k')
  plot(x_hist(:, 1), x_hist(:, 4), '.--r')
  axis([0 1 0 1]*1000/2)
  grid on


end
