%% IMM - Trabalho final

curva_suave = load('curva_suave');
manobra = load('manobra');

%% Estado inicial e modelo
dt = 4;                                        % Intervalo de amostragem
x0 = [3; 0; 0]*1e4;                           % Média inicial (escolhemos -3e4 para
                                               % posição inicial e zero para velocidade inicial)
P0 = 1e7*eye(length(x0));                      % Matriz de covariância inicial (Alta Incerteza)

  F1 = [1 dt 0;                                     % Modelo de posicao quase constante
      0 1  0;
      0 0  0];
F2 = [1 dt dt.^2/2; ...                                % Modelo de velocidade quase constante
      0 1  dt;
      0 0  1];
F3 = [1 0 0;                                     % Modelo de posicao quase constante
      0 0  0;
      0 0  0];

Q1 = [0 0 0; ...                                  % Modelo de ruído somente na componente de velocidade
      0 1 0; ...
      0 0 0];
Q2 = [0 0 0; ...                              % Modelo de ruído somente na componente de velocidade
      0 0 0;
      0 0 1]
Q3 = [1 0 0; ...                                  % Modelo de ruído somente na componente de velocidade
      0 0 0; ...
      0 0 0];

H1 = [1 0 0] ;                              % Medimos somente a posição (Estados 1 and 3)
H2 = [1 0 0];  ...                             % Medimos somente a posição (Estados 1 and 3)
H3 = [1 0 0] ;                              % Medimos somente a posição (Estados 1 and 3)

sigma_r = sqrt(5e8);
R = sigma_r^2*eye(size(H1, 1));                 % Matriz de covariância do ruído

M = [0.90 0.05 0.05; ...
     0.05 0.90 0.05;
     0.05 0.05 0.90];
mu = [0.5 0.5 0.5];
zv = Xevol(1,:)' + sigma_r*randn(size(Xevol(1, :)'));

Q1 = Q1*1e4;
Q2 = Q2*1e4;
Q3 = Q3*1e4;

[x_est, mu_est] = imm3(x0, P0, zv, F1, H1, Q1, F2, H2, Q2, F3, H3, Q3, R, mu, M);

%% Plotting curves
figure(1)
plot(zv, 'o')
hold on
plot(x_est(1, :), '.--')
grid

figure(2)
plot(mu_est', '.--')
grid
