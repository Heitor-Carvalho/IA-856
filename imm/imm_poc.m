%% Trabalho final - Interacting Multiple Model

addpath('../kalman')

reta       = load('../data_ex/reta');
curva      = load('../data_ex/curva_suave');
manobra    = load('../data_ex/manobra');
zv_reta    = reta.x_m';
zv_curva   = curva.x_m';
zv_manobra = manobra.x_m';

%% Estado inicial e modelo
dt = 4;                                        % Intervalo de amostragem
x0 = [-3; 0; 0; -3; 0; 0]*1e4;                 % Média inicial (escolhemos -3e4 para
                                               % posição inicial e zero para velocidade inicial)
P0 = 1e6*eye(length(x0));                      % Matriz de covariância inicial (Alta Incerteza)

F1 = [1 dt 0;                                  % Modelo de velocidade quase constante
      0 1  0;
      0 0  0];
F2 = [1 dt dt.^2/2; ...                        % Modelo de aceleração quase constante
      0 1  dt;
      0 0  1];

Q1 = [0 0 0 0 0 0; ...                         % Modelo de ruído somente na componente de velocidade
      0 1 0 0 0 0; ...
      0 0 0 0 0 0; ...
      0 0 0 0 0 0; ...
      0 0 0 0 1 0; ...
      0 0 0 0 0 0];
Q2 = [0 0 0 0 0 0; ...                         % Modelo de ruído somente na componente de aceleração
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

%% Comparando o IMM com o filtro de Kalman (caso reta)
M1 = [0.995 0.005; 0.005 0.995];
M2 = [0.98 0.02; 0.02 0.98];
M3 = [0.95 0.05; 0.05 0.95];

[x_est_reta_m1, ~] = imm2(x0, P0, zv_reta, FF1, H1, Q1, FF2, H2, Q2, R, mu, M1);
[x_est_reta_m2, ~] = imm2(x0, P0, zv_reta, FF1, H1, Q1, FF2, H2, Q2, R, mu, M2);
[x_est_reta_m3, ~] = imm2(x0, P0, zv_reta, FF1, H1, Q1, FF2, H2, Q2, R, mu, M2);

[x_kalman_reta.x_est, ~, ~, ~, ~] = kalman_filter_wrapper(x0, P0, zv_reta, FF1, H1, Q1, R);

start_sample = 15;
mse_imm_m1 = mean( abs(x_est_reta_m1(1, :)-reta.x_c(:, 1)').^2 + abs(x_est_reta_m1(4, :) - reta.x_c(:, 2)').^2 );
mse_imm_m2 = mean( abs(x_est_reta_m2(1, :)-reta.x_c(:, 1)').^2 + abs(x_est_reta_m2(4, :) - reta.x_c(:, 2)').^2 );
mse_imm_m3 = mean( abs(x_est_reta_m3(1, :)-reta.x_c(:, 1)').^2 + abs(x_est_reta_m3(4, :) - reta.x_c(:, 2)').^2 );
mse_kalman = mean( abs(x_kalman_reta.x_est(1, :)-reta.x_c(:, 1)').^2 + abs(x_kalman_reta.x_est(4, :) - reta.x_c(:, 2)').^2 );

mse_str = sprintf('|  %.2f  |  %.2f  |  %.2f  |  %.2f  |', mse_kalman, mse_imm_m1, mse_imm_m2, mse_imm_m3)

disp('|   Kalman     |  IMM - M1   |   IMM - M2  |  IMM - M3 |')
disp(mse_str)

% Como era esperado para este caso, o MSE do algorítmo IMM ficou sempre
% acima do obtido pelo filtro de kalman, sendo este o limite inferior 
% para este caso.
% Também é possível ver que neste caso, ao usar matrizes cuja probabilidade de 
% transição entre os modelos é mais baixa, o MSE obtido pelo algorítmo IMM fica 
% mais próximo ao obtido pelo filtro de kalman. 
% Isto se deve ao fato do filtro se comportar mais como o modelo 1
% (velocidade constante) ao diminuirmos a probabilidade de transição de
% estado.

%% Comentários e conclusão
% Neste trabalho foi possível observar que ao usar múltiplos modelos com 
% seus parâmetros ajustados adequadamente, pudemos alcançar um ótimo
% trade-off entre filtragem e tracking. Verificamos que os pesos dados a
% cada modelo para estimativa final (correspondentes a probabilidade de
% cada modelo), foram coerentes com o esperado. No caso onde o objeto se
% movimentou com velocidade constante, pudemos ver que o modelo de
% velocidade quase constante teve maior probabilidade. Já na situação de
% manobra pode-se ver que as probabilidades se inverteram e o modelo de
% aceleração quase constante teve maior peso na estimativa final. 

%%
close all