function [x_est, mu_est] = imm2(x0, P0, zv, F1, H1, Q1, F2, H2, Q2, R, mu, M)

  x1 = x0;                                     % Loading initial state mean 1
  P1 = P0;                                     % Loading initial cov matrix 1
  x2 = x0;                                     % Loading initial state mean 2
  P2 = P0;                                     % Loading initial cov matrix 2

  x_est = zeros(size(x0,1), size(zv, 1));
  mu_est = zeros(size(mu,2), size(zv, 1));

  for i=1:length(zv)
    z = zv(i, :)';

    % 1 - Calcula das probabilidades de mistura
    % Neste etapa calculamentos a probabilidade das 
    % transições entre estados dado a medida - u(k-1/k-1)

    cbar = mu*M;
    % cbar é probabilidade total de que o sistema esteja no 
    % no modelo j (Estado da cadeia de Markov). 
    % Estes valores serviram com a constante de normalização
    

    omega = zeros(size(M));
    for j=1:length(M)
      for k=1:length(M)
        omega(j,k) = M(j, k)*mu(j)/cbar(k);
      end
    end

    % 2. Mistura

    % Misturamos as estimativas com base na probabilidade 
    % de transição de estado
    w1 = omega(:, 1);
    x1mix = zeros(length(x0), 1);
    x1mix = x1*w1(1) + x2*w1(2);
    P1mix = zeros(size(P0));
    y1 = x1-x1mix;

    w2 = omega(:, 2);
    x2mix = zeros(length(x0), 1);
    x2mix = x1*w2(1) + x2*w2(2);
    P2mix = zeros(size(P0));
    y2 = x2-x2mix;

    P1mix = w1(1)*(y1*y1' + P1) + w1(2)*(y2*y2' + P2);
    P2mix = w2(1)*(y1*y1' + P1) + w2(2)*(y2*y2' + P2);

    % 3 Mode-matched filtering

    % Passo de predição do filtro 1 - Com as estimativas misturadas
    [x1, P1] = kalman_predict(F1, Q1, 0, 0, x1mix, P1mix);

    % Passo de atualização do filtro 1 - Com as estimativas misturadas
    [x1, P1,residual1, S1, K] = kalman_update(H1, R, z, x1, P1);

    % Passo de predição do filtro 2 - Com as estimativas misturadas
    [x2, P2] = kalman_predict(F2, Q2, 0, 0, x2mix, P2mix);

    % Passo de atualização do filtro 2 - Com as estimativas misturadas
    [x2, P2, residual2, S2, K] = kalman_update(H2, R, z, x2, P2);

    % Likelihood
    L =  [mvnpdf(zeros(size(residual1)), residual1, S1), mvnpdf(zeros(size(residual2)), residual2, S2)];


    % 4 Atualização da probabilidade de cada modo
    % Equação de Bayes - Prio*Likelihood/Normalization
    mu = cbar.*L;
    mu = mu/sum(mu, 2);

    % 5 Estimativa dos estados e covariância.

    x = zeros(length(x0), 1);
    P = zeros(2,2);
    x = x1*mu(1) + x2*mu(2);
    y1 = x1-x;
    y2 = x2-x;
    P = mu(1)*(y1*y1' + P1) + mu(2)*(y2*y2' + P2);


    x_est(:, i) = x;
    mu_est(:, i) = mu;
  end

end

