%% IMM - Trabalho final

reta = load('data_ex/reta');
curva = load('data_ex/curva_suave');
manobra = load('data_ex/manobra');
zv_reta = reta.x_m';
zv_curva = curva.x_m';
zv_manobra = manobra.x_m';

%% Estado inicial e modelo
dt = 4;                                        % Intervalo de amostragem
x0 = [-3; 0; 0; -3; 0; 0]*1e4;                 % Média inicial (escolhemos -3e4 para
                                               % posição inicial e zero para velocidade inicial)
P0 = 1e4*eye(length(x0));                      % Matriz de covariância inicial (Alta Incerteza)

F1 = [1 dt 0;                                  % Modelo de velocidade quase constante
      0 1  0;
      0 0  0];
F2 = [1 dt dt.^2/2; ...                        % Modelo de posição quase constante
      0 1  dt;
      0 0  1];

Q1 = [0 0 0 0 0 0; ...                         % Modelo de ruído somente na componente de velocidade
      0 1 0 0 0 0; ...
      0 0 0 0 0 0; ...
      0 0 0 0 0 0; ...
      0 0 0 0 1 0; ...
      0 0 0 0 0 0];
Q2 = [0 0 0 0 0 0; ...                         % Modelo de ruído somente na componente de velocidade
      0 0 0 0 0 0; ...
      0 0 1 0 0 0; ...
      0 0 0 0 0 0; ...
      0 0 0 0 0 0; ...
      0 0 0 0 0 1];

H1 = [1 0 0 0 0 0; ...                         % Medimos somente a posição (Estados 1 and 3)
      0 0 0 1 0 0];
H2 = [1 0 0 0 0 0;  ...                        % Medimos somente a posição (Estados 1 and 3)
      0 0 0 1 0 0];

sigma_r = 1200;
R = sigma_r^2*eye(size(H1, 1));                % Matriz de covariância do ruído

FF1 = [F1 zeros(3); zeros(3) F1];              % Modelo expandido com dois modelos de velocidade
                                               % quase constante independentes resultando em
                                               % uma matriz 4x4

FF2 = [F2 zeros(3); zeros(3) F2];              % Modelo expandido com dois modelos de velocidade
                                               % quase constante independentes resultando em
                                               % uma matriz 4x4

M = [0.98 0.02; ...
     0.02 0.98];
mu = [0.5 0.5];

Q1 = Q1*0.05;
Q2 = Q2*0.08;

%% Filtering trajectories with IMM
[x_est_reta, mu_est_reta] = imm2(x0, P0, zv_reta, FF1, H1, Q1, FF2, H2, Q2, R, mu, M);
[x_est_curva, mu_est_curva] = imm2(x0, P0, zv_curva, FF1, H1, Q1, FF2, H2, Q2, R, mu, M);
[x_est_manobra, mu_est_manobra] = imm2(x0, P0, zv_manobra, FF1, H1, Q1, FF2, H2, Q2, R, mu, M);

%% Plotting curves - Reta
figure(1)
plot(reta.x_m(1, :), reta.x_m(2, :), 'o')
hold on
plot(x_est_reta(1, :), x_est_reta(4, :), '.--')
legend('Medidas', 'Sinal Filtrado', 'Location', 'SouthEast')
grid

figure(2)
plot(mu_est_reta', '.--')
legend('Velocidade Constante', 'Aceleração Constante', 'Location', 'SouthOutside')
grid

%% Plotting curves - Curva suave
figure(3)
plot(curva.x_m(1, :), curva.x_m(2, :), 'o')
hold on
plot(x_est_curva(1, :), x_est_curva(4, :), '.--')
legend('Medidas', 'Sinal Filtrado', 'Location', 'SouthEast', 'Location', 'SouthEast')
grid

figure(4)
plot(mu_est_curva', '.--')
legend('Velocidade Constante', 'Aceleração Constante', 'Location', 'SouthOutside')
grid

%% Plotting curves - Manobra

figure(5)
plot(manobra.x_m(1, :), manobra.x_m(2, :), 'o')
hold on
plot(x_est_manobra(1, :), x_est_manobra(4, :), '.--')
legend('Medidas', 'Sinal Filtrado', 'Location', 'SouthEast', 'Location', 'SouthEast')
grid

figure(6)
plot(mu_est_manobra', '.--')
legend('Velocidade Constante', 'Aceleração Constante', 'Location', 'SouthOutside')
grid

%%
close all