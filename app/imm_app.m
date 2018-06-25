%% IMM Kalman Filter App
clear all
close all

%% Declaring global variable
global F1
global H1
global Q1
global F2
global H2
global Q2
global pos_hist
global pos_noise_hist
global x_hist
global mu_hist
global sigma_r
global R
global x1
global mu
global M
global P1
global x2
global P2

%% Estado inicial e modelo
dt = 0.8;                                          % Intervalo de amostragem
x0 = [0; 0; 0; 0; 0; 0]*1e3;                       % Média inicial (escolhemos -3e4 para
                                                   % posição inicial e zero para velocidade inicial)
P0 = 1e8*eye(length(x0));                          % Matriz de covariância inicial (Alta Incerteza)

F1 = [1 dt 0;                                      % Modelo de velocidade quase constante
      0 1  0;
      0 0  0];
F2 = [1 dt dt.^2/2; ...                            % Modelo de aceleração quase constante
      0 1  dt;
      0 0  1];

Q1 = [0 0 0 0 0 0; ...                             % Modelo de ruído somente na componente de velocidade
      0 1 0 0 0 0; ...
      0 0 0 0 0 0; ...
      0 0 0 0 0 0; ...
      0 0 0 0 1 0; ...
      0 0 0 0 0 0];
Q2 = [0 0 0 0 0 0; ...                             % Modelo de ruído somente na componente de aceleração
      0 0 0 0 0 0; ...
      0 0 1 0 0 0; ...
      0 0 0 0 0 0; ...
      0 0 0 0 0 0; ...
      0 0 0 0 0 1];

H1 = [1 0 0 0 0 0; ...                              % Medimos somente a posição (Estados 1 and 3)
      0 0 0 1 0 0];
H2 = [1 0 0 0 0 0;  ...                             % Medimos somente a posição (Estados 1 and 3)
      0 0 0 1 0 0];


F1 = [F1 zeros(3); zeros(3) F1];                    % Modelo expandido com dois modelos de velocidade
                                                    % quase constante independentes resultando em
                                                    % uma matriz 4x4

F2 = [F2 zeros(3); zeros(3) F2];                    % Modelo expandido com dois modelos de aceleração
                                                    % quase constante independentes resultando em
                                                    % uma matriz 4x4
sigma_r = 8;
R = sigma_r^2*eye(size(H1, 1));                     % Matriz de covariância do ruído

sigma_q1 = 0.1;
sigma_q2 = 5;
Q1 = sigma_q1*Q1;
Q2 = sigma_q2*Q2;

M = [0.98 0.02; ...
     0.02 0.98];

Q1 = Q1*0.5;
Q2 = Q2*0.8;

mu = [0.5 0.5];

hist_len       = 25;
pos_hist       = zeros(hist_len, 2);
pos_noise_hist = zeros(hist_len, 2);
x_hist         = zeros(hist_len, 6);
mu_hist        = zeros(hist_len, 2);
x1 = x0;
P1 = P0;
x2 = x0;
P2 = P0;

t = timer('StartDelay', 0, 'Period', 0.1, 'ExecutionMode', 'fixedRate');
t.TimerFcn = @ImmMoveCallback;
start(t)
