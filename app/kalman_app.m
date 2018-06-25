%% Kalman Filter App
clear all
close all

%% Declaring global variable
global F
global H
global Q
global R
global pos_hist
global pos_noise_hist
global x_hist
global lk_hist
global sigma_r
global x
global P

%% Estado inicial e modelo
dt = 1;                                        % Intervalo de amostragem
x0 = [0; 0; 0; 0]*1e3;                       % Média inicial (escolhemos -3e4 para
                                               % posição inicial e zero para velocidade inicial)
P0 = 1e9*eye(length(x0));                      % Matriz de covariância inicial (Alta Incerteza)

F = [1 dt; ...                                 % Modelo de velocidade quase constante
     0 1];
Q = [0 0 0 0; ...                              % Modelo de ruído somente na componente de velocidade
     0 1 0 0; ...
     0 0 0 0; ...
     0 0 0 1];

H = [1 0 0 0;  ...                             % Medimos somente a posição (Estados 1 and 3)
     0 0 1 0];

sigma_r = 8;
R = sigma_r^2*eye(size(H, 1));                 % Matriz de covariância do ruído

sigma_q = 0.5;
Q = sigma_q*Q;

F = [F zeros(2); zeros(2) F];                 % Modelo expandido com dois modelos de velocidade
                                              % quase constante independentes resultando em
                                              % uma matriz 4x4
hist_len       = 25;
pos_hist       = zeros(hist_len, 2);
pos_noise_hist = zeros(hist_len, 2);
x_hist         = zeros(hist_len, 4);
lk_hist        = zeros(hist_len, 1);
x = x0;
P = P0;

t = timer('StartDelay', 0, 'Period', 0.1, 'ExecutionMode', 'fixedRate');
t.TimerFcn = @MoveCallback;
start(t)