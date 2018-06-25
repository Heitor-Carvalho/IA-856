%% Lista - Filtro de Kalman
clear all

%% Carregando dados da trajetória
reta        = load('./data_ex/reta');
curva_suave = load('./data_ex/curva_suave');
manobra     = load('./data_ex/manobra');

%% Estado inicial e modelo
dt = 4;                                        % Intervalo de amostragem
x0 = [-3; 0; -3; 0]*1e4;                       % Média inicial (escolhemos -3e4 para
                                               % posição inicial e zero para velocidade inicial)
P0 = 1e4*eye(length(x0));                      % Matriz de covariância inicial (Alta Incerteza)

F = [1 dt; ...                                 % Modelo de velocidade quase constante
     0 1];
Q = [0 0 0 0; ...                              % Modelo de ruído somente na componente de velocidade
     0 1 0 0; ...
     0 0 0 0; ...
     0 0 0 1];

H = [1 0 0 0;  ...                             % Medimos somente a posição (Estados 1 and 3)
     0 0 1 0];

sigma_r = 1200;
R = sigma_r^2*eye(size(H, 1));                 % Matriz de covariância do ruído


FF = [F zeros(2); zeros(2) F];                 % Modelo expandido com dois modelos de velocidade
                                               % quase constante independentes resultando em
                                               % uma matriz 4x4

%% Variando o rúido do sistema (matriz Q)

q_sigma = [1 10 50 100];
q_sig_data_reta = cell(length(q_sigma), 1);
q_sig_data_curva_suave = cell(length(q_sigma), 1);
q_sig_data_manobra = cell(length(q_sigma), 1);
for i =1:length(q_sigma)
  Qt = Q*q_sigma(i).^2;
  [q_sig_data_reta{i}.x_est        , q_sig_data_reta{i}.likelyhood, ...
   q_sig_data_reta{i}.residual_hist, q_sig_data_reta{i}.trace_hist, q_sig_data_reta{i}.P] = kalman_filter_wrapper(x0, P0, reta.x_m', FF, H, Qt, R);
  [q_sig_data_curva_suave{i}.x_est        , q_sig_data_curva_suave{i}.likelyhood, ...
   q_sig_data_curva_suave{i}.residual_hist, q_sig_data_curva_suave{i}.trace_hist, q_sig_data_curva_suave{i}.P] = kalman_filter_wrapper(x0, P0, curva_suave.x_m', FF, H, Qt, R);
  [q_sig_data_manobra{i}.x_est        , q_sig_data_manobra{i}.likelyhood, ...
   q_sig_data_manobra{i}.residual_hist, q_sig_data_manobra{i}.trace_hist, q_sig_data_manobra{i}.P] = kalman_filter_wrapper(x0, P0, manobra.x_m', FF, H, Qt, R);
end

%% Trajetória estimada para diferentes valores de $\sigma_q^2$

figure(1)
plot(reta.x_m(1, :), reta.x_m(2, :), 'o')
hold on
for i = 1:length(q_sigma)
  plot(q_sig_data_reta{i}.x_est(1, :)', q_sig_data_reta{i}.x_est(3, :)', '.--')
end
legend('Noisy Measurments', 'Estimation - \sigma_q = 1', 'Estimation - \sigma_q = 10', ...
       'Estimation - \sigma_q = 50', 'Estimation - \sigma_q = 100', 'Location', 'SouthEast')
xlabel('Estado 1')
ylabel('Estado 3')
grid

figure(2)
plot(curva_suave.x_m(1, :), curva_suave.x_m(2, :), 'o')
hold on
for i = 1:length(q_sigma)
  plot(q_sig_data_curva_suave{i}.x_est(1, :)', q_sig_data_curva_suave{i}.x_est(3, :)', '.--')
end
legend('Noisy Measurments', 'Estimation - \sigma_q = 1', 'Estimation - \sigma_q = 10', ...
       'Estimation - \sigma_q = 50', 'Estimation - \sigma_q = 100', 'Location', 'SouthEast')
xlabel('Estado 1')
ylabel('Estado 3')
grid

figure(3)
plot(manobra.x_m(1, :), manobra.x_m(2, :), 'o')
hold on
for i = 1:length(q_sigma)
  plot(q_sig_data_manobra{i}.x_est(1, :)', q_sig_data_manobra{i}.x_est(3, :)', '.--')
end
legend('Noisy Measurments', 'Estimation - \sigma_q = 1', 'Estimation - \sigma_q = 10', ...
       'Estimation - \sigma_q = 50', 'Estimation - \sigma_q = 100', 'Location', 'SouthEast')
xlabel('Estado 1')
ylabel('Estado 3')
grid

%% Trajetória estimada para diferentes valores de $\sigma_q^2$ - Zoon-in
figure(4)
plot(reta.x_m(1, :), reta.x_m(2, :), 'o')
hold on
for i = 1:length(q_sigma)
  plot(q_sig_data_reta{i}.x_est(1, :)', q_sig_data_reta{i}.x_est(3, :)', '.--', 'LineWidth', 2.5)
end
legend('Noisy Measurments', 'Estimation - \sigma_q = 1', 'Estimation - \sigma_q = 10', ...
       'Estimation - \sigma_q = 50', 'Estimation - \sigma_q = 100', 'Location', 'SouthEast')
xlabel('Estado 1')
ylabel('Estado 3')
xlim([0 16000])
grid

figure(5)
plot(curva_suave.x_m(1, :), curva_suave.x_m(2, :), 'o')
hold on
for i = 1:length(q_sigma)
  plot(q_sig_data_curva_suave{i}.x_est(1, :)', q_sig_data_curva_suave{i}.x_est(3, :)', '.--', 'LineWidth', 2.5)
end
legend('Noisy Measurments', 'Estimation - \sigma_q = 1', 'Estimation - \sigma_q = 10', ...
       'Estimation - \sigma_q = 50', 'Estimation - \sigma_q = 100', 'Location', 'SouthEast')
xlabel('Estado 1')
ylabel('Estado 3')
xlim([-12000 -1000])
grid

figure(6)
plot(manobra.x_m(1, :), manobra.x_m(2, :), 'o')
hold on
for i = 1:length(q_sigma)
  plot(q_sig_data_manobra{i}.x_est(1, :)', q_sig_data_manobra{i}.x_est(3, :)', '.--', 'LineWidth', 2.5)
end
legend('Noisy Measurments', 'Estimation - \sigma_q = 1', 'Estimation - \sigma_q = 10', ...
       'Estimation - \sigma_q = 50', 'Estimation - \sigma_q = 100', 'Location', 'SouthEast')
xlabel('Estado 1')
ylabel('Estado 3')
xlim([-10000 4000])
grid

% * Como pode ser visto na figura para o caso reta, a melhor trajetória foi
% obtida quando $\sigma_q$ = 1, ou seja, com o menor valor de $\sigma_q$ testado. Este resultado é
% esperado pois neste caso um modelo de velocidade constante (reta) descreve de maneira correta nosso
% sistema e ao fazer $\sigma_q$ pequeno, estamos nos aproximando mais deste modelo.
% * Já para os casos de trajetória com curva suave e manobra, pode-se notar
% que $\sigma_q$ = 1 faz com que o filtro não altere sua estimativa de maneira
% rápida o suficiente para acompanhar a trajetória. Para que o filtro acompanhe a trajetória
% corretamente, foi necessário aumentar $\sigma_q$ deixando o modelo mais próximo de
% um modelo de velocidade variável. Após alguns testes chegou-se ao valor
% de $\sigma_q$ igual 50 para um valor que funcionasse tanto para a trajetória de
% curva suave quanto para manobra.
% * Também se pode notar que para $\sigma_q$ muito elevado, nosso modelo possui
% muita flexibilidade e com isso acaba seguindo o rúido, isto é, o valor estimado pelo
% filtro é pratimente o valor medido. Este fato também pode ser visto nas figuras
% observando os valores estimados para $\sigma_q$ elevado.
% * Com isso, vemos que nosso filtro não poderá ter o melhor valor de $\sigma_q$
% para as três trajetórias propostas e que se quisermos um filtro que
% funcione nas três situações iremos precisar escolher o menor valor de $\sigma_q$
% que resulte em um resultado aceitável para a trajetória de manobra, pois esta exige um
% $\sigma_q$ mais elevado.
% * Aqui, escolhemos o valor de $\sigma_q$ = 50 como valor para usarmos nas três trajetórias,
% porém pode-se notar que este valor não irá produzir o melhor resultado quando nosso
% alvo se more em uma linha reta (O que pode acontecer boa parte do tempo), ou fizer uma curva mais suave.

%% MSE dos estados de trabalho

mse_q_reta        = zeros(1, length(q_sigma));
mse_q_curva_suave = zeros(1, length(q_sigma));
mse_q_manobra     = zeros(1, length(q_sigma));
for i = 1:length(q_sigma)
  mse_q_reta(i)        = mean( (reta.x_c(:,1) - q_sig_data_reta{i}.x_est(1, :)').^2 + (reta.x_c(:,2) - q_sig_data_reta{i}.x_est(3, :)').^2 );
  mse_q_curva_suave(i) = mean( (curva_suave.x_c(:,1) - q_sig_data_curva_suave{i}.x_est(1, :)').^2 + (curva_suave.x_c(:,2) - q_sig_data_curva_suave{i}.x_est(3, :)').^2 );
  mse_q_manobra(i)     = mean( (manobra.x_c(:,1) - q_sig_data_manobra{i}.x_est(1, :)').^2 + (manobra.x_c(:,2) - q_sig_data_manobra{i}.x_est(3, :)').^2 );
end

disp('Information - Caso reta')
for i = 1:length(q_sigma)
  info = sprintf('MSE = %.3g | P(1,1) = %.3g | P(3,3) = %.3g | Q = %.3g | R = %.3g', mse_q_reta(i), q_sig_data_reta{i}.P(1,1), ...
                                                                                     q_sig_data_reta{i}.P(3,3), q_sigma(i), sigma_r.^2);
  disp(info)
end

disp('Information - Caso curva_suave')
for i = 1:length(q_sigma)
  info = sprintf('MSE = %.3g | P(1,1) = %.3g | P(3,3) = %.3g | Q = %.3g | R = %.3g', mse_q_reta(i), q_sig_data_curva_suave{i}.P(1,1), ...
                                                                                     q_sig_data_curva_suave{i}.P(3,3), q_sigma(i), sigma_r.^2);
  disp(info)
end

disp('Information - Caso manobra')
for i = 1:length(q_sigma)
  info = sprintf('MSE = %.3g | P(1,1) = %.3g | P(3,3) = %.3g | Q = %.3g | R = %.3g', mse_q_reta(i), q_sig_data_manobra{i}.P(1,1), ...
                                                                                     q_sig_data_manobra{i}.P(3,3), q_sigma(i), sigma_r.^2);
  disp(info)
end

% Aqui temos um resumo das informações dos filtros avaliados onde podemos
% verificar a ordem de grandez as variáveis Q, R, das variâncias dos
% estados 1 e 3 e do MSE calculado em relação aos estados reais.
% Pode-se ver que os valores de variância estimado pelo filtro (matriz P)
% foram bem próximos aos valores de MSE em relação ao estado real. Isso
% mostra que nosso filtro possui uma boa capacidade de estimar não somente
% a média, mas também a incerteza em relação ao valor estimado. É
% importante notar que essa precisão só é válida caso as hipóteses de nosso
% filtro sejam satisfeitas (ou próximas).
% Também notamos que o aumento de Q fez com que o MSE aumentasse se
% aproximando do valor de R. Isso é esperado, pois ao aumentar Q fazemos
% com que nossas estimativas dos estados siguam os valores medidos, dessa
% maneira nosso filtro terá um MSE mais próximo a variância das medidas, ou
% seja, igual a variância de R.

%% Variando o estado inicial e a matriz de covariância
% Para mostrar o efeito de variação do estado inicial e
% da covariância inicial, usaremos somente a trajetória de curva suave.
% Esta decisão foi escolhida pois este caso já é suficiente para demonstrar
% o comportamento esperado e assim iremos evitar análises repetidas.
x0 = [-2; 0; 1; 0]*1e4;                       % Posicição inicial
Qt = Q*q_sigma(3).^2;

sigma_p0 = [0.001, 100, 1000];

p_sig_data_curva_suave = cell(length(q_sigma), 1);
for i =1:length(sigma_p0)
  Pt = sigma_p0(i).^2*eye(length(x0)); 
  [p_sig_data_curva_suave{i}.x_est        , p_sig_data_curva_suave{i}.likelyhood, ...
   p_sig_data_curva_suave{i}.residual_hist, p_sig_data_curva_suave{i}.trace_hist, p_sig_data_curva_suave{i}.P] = kalman_filter_wrapper(x0, Pt, curva_suave.x_m', FF, H, Qt, R);
end


figure(7)
plot(curva_suave.x_m(1, :), curva_suave.x_m(2, :), 'o')
hold on
for i = 1:length(sigma_p0)
  plot(p_sig_data_curva_suave{i}.x_est(1, :)', p_sig_data_curva_suave{i}.x_est(3, :)', '^--')
end
legend('Noisy Measurments', 'Estimation - \sigma_{x0} = 0.001', 'Estimation - \sigma_{x0} = 100', ...
       'Estimation - \sigma_{x0} = 1000', 'Location', 'NorthWest')
xlabel('Estado 1')
ylabel('Estado 3')
xlim([-3.5 -1]*1e4)
grid

% Aqui, podemos ver que ao iniciar a matriz de covariância do estado com um
% valor pequeno fazemos com que seja necessário mais passos até que estado
% passse a acompanhar a trajetória dada pelos valores medidos. Isto faz
% sentido pois estamos informando ao nosso filtro que há pouca incerteza
% sobre o estado do sistema e que nosso chute inicial é bom. Já ao aumentar
% a variância da matriz de covariância indicamos uma alta incerteza sobre a
% posição inicial, o que faz com que o estado inicial passe a acompanhar as
% medidas mais rapidamente. Novamente esse comportamente é esperado, pois estamos
% informando o filtro de que nossa inicialização possui muita incerteza e
% assim o filtro irá rapidamente abandonar nossa estimativa inicial e
% acompanhar as medidas.

    
%% Estabilidade do filtro
% Vamos demonstrar a estabilidade do filtro considerando o caso em que 
% $\sigma_q$ = 50 e a trajotória é uma reta, para isso iremos analisar o 
% valor da matriz P e comparar com valor obtido através da solução da 
% equação algébrica de Ricatti. Vamos também provar as condições necessárias 
% e suficiente para a estabilidade do filtro

% Verificando se o par (F,H) é observável, neste caso a matriz de
% observabilidade deve ter rank = 4
observability_matrix = obsv(FF, H);
obsv_rank = rank(observability_matrix)

% Verificando se o par (F,Q^0.5) é controlável, neste caso a matriz de
% controlabilidade deve ter rank = 4
controlability_matrix = ctrb(FF, Qt^0.5);
ctrb_rank = rank(controlability_matrix)

% Com base nas condições acima e que P0 > 0, temos que o filtro é estável e
% tem solução unica. Podemos então calcular a solução estacionária com o seguinte
% comando e comparar com a matriz de covariância do filtro P
P_stationary = dare(FF', H', Qt, R)           % Matriz de covariância estacionária, o comando dare
                                              % resolve a equação de algébrica de Ricatti
P = q_sig_data_reta{3}.P                      % Matriz de covaruância após o passo de predição

%% Análise da verossimilhança na estimação
% Para mostrar a variação da verossimilhança, vamos usar uma
% valor de $\sigma_q$ = 10 e observar a verossimilhança para os
% trajetórias retas e manobra.

figure(8)
subplot(2, 1, 1)
plot(q_sig_data_reta{2}.likelyhood)
ylim([0 9e-8])
legend('Trajetória - Reta', 'Location', 'SouthWest')
ylabel('Verossimilhança')
xlabel('Samples')
grid
subplot(2, 1, 2)
plot(q_sig_data_manobra{2}.likelyhood)
ylim([0 9e-8])
legend('Trajetória - Manobra', 'Location', 'SouthWest')
ylabel('Verossimilhança')
xlabel('Samples')
grid


% Antes de comentar o comportamento da verossimilhança é preciso 
% lembrar que esta é calculada usando o resíduo entre a estimativa do modelo, 
% a medida processada e também a matriz de covarância das medidas (matriz
% S). Isto explica o comportamento oscilátório da verossimilhança já que
% temos ruídos em nossa estimativa.
% Ao analisar as curvas da verossimilhança para as trajetórias de reta e manobra 
% para um valor de $sigma_q$ baixo, notamos que para o caso de linha reta os valores 
% apesar de ruídosos sempre se mantiveram próximos a um valor médio. Porém, para a 
% trajetória em manobra é possível notar uma grande queda no valor da verssimilhança
% durante as amostras 45 e 60, esta queda corresponde ao período onde o filtro não
% consegue acompanhar os valores das medidas devido ao baixo $sigma_q$ usado. 
% Ou seja, podemos usar a verossimilhança das medidas para estimar se nosso filtro 
% está funcionando adequadamente, não somente isso, mas também podemos usar essa 
% informação para ativar uma alteração ou adaptação de parâmetros de nosso filtro. 

%% Proposta de adaptação do ruído Q
% Com base nessa análise, partimos para o último passo na análise do filtro
% de kalman, vamos tentar melhorar nossa filtragem usando um $\sigma_q$
% pequeno e o aumentando quando nossa verossimilhanção for menor 1e-8

x0 = [-3; 0; -3; 0]*1e4;                       % Média inicial (escolhemos -3e4 para
                                               % posição inicial e zero para velocidade inicial)
P0 = 1e4*eye(length(x0));                      % Matriz de covariância inicial (Alta Incerteza)
Qmin = Q*1.^2;                                 % $sigma_q$ inicial
Qmax = Q*50.^2;                                % $sigma_q$ para casos onde verossimilhança menor que 1e-8
limit = 1e-8;
[q_sig_data_reta_adpt.x_est         , q_sig_data_reta_adpt.likelyhood, ...
  q_sig_data_reta_adpt.residual_hist, q_sig_data_reta_adpt.trace_hist, q_sig_data_reta_adpt.P] =  ...
  kalman_filter_wrapper_adpt(x0, P0, reta.x_m', FF, H, Qmin, Qmax, R, limit);
[q_sig_data_curva_suave_adpt.x_est         , q_sig_data_curva_suave_adpt.likelyhood, ...
  q_sig_data_curva_suave_adpt.residual_hist, q_sig_data_curva_suave_adpt.trace_hist, q_sig_data_curva_suave_adpt.P] =  ...
  kalman_filter_wrapper_adpt(x0, P0, curva_suave.x_m', FF, H, Qmin, Qmax, R, limit);
[q_sig_data_manobra_adpt.x_est         , q_sig_data_manobra_adpt.likelyhood, ...
  q_sig_data_manobra_adpt.residual_hist, q_sig_data_manobra_adpt.trace_hist, q_sig_data_manobra_adpt.P] =  ...
  kalman_filter_wrapper_adpt(x0, P0, manobra.x_m', FF, H, Qmin, Qmax, R, limit);

figure(9)
plot(reta.x_m(1, :), reta.x_m(2, :), 'o')
hold on 
plot(q_sig_data_reta_adpt.x_est(1, :), q_sig_data_reta_adpt.x_est(3, :), '.--')
xlabel('Estado 1')
ylabel('Estado 3')
grid

figure(10)
plot(curva_suave.x_m(1, :), curva_suave.x_m(2, :), 'o')
hold on 
plot(q_sig_data_curva_suave_adpt.x_est(1, :), q_sig_data_curva_suave_adpt.x_est(3, :), '.--')
xlabel('Estado 1')
ylabel('Estado 3')
grid

figure(11)
plot(manobra.x_m(1, :), manobra.x_m(2, :), 'o')
hold on 
plot(q_sig_data_manobra_adpt.x_est(1, :), q_sig_data_manobra_adpt.x_est(3, :), '.--')
xlabel('Estado 1')
ylabel('Estado 3')
grid

% Como podemos ver pelas figuras, a alteração de $sigma_q$ de 1 para 50
% quando a verossimilhança cai para baixo de 1e-8 fez com que nossa
% estimativa fosse mais suave e próxima ao valor esperado em todas três
% trajetórias. Assim, temos um forma de combater o trade-off enfrentado na
% hora da escolha do valor de Q. Por fim, vamos mostrar o comportamento da
% verossimilhança para a trajetória em manóbra e ver analisar seu
% comportamento.

figure(12)
plot(q_sig_data_manobra_adpt.likelyhood, '.--')
legend('Trajetória - Manobra', 'Location', 'SouthWest')
ylabel('Verossimilhança')
xlabel('Samples')
grid

% Podemos ver que diferente do caso sem adaptação de Q, a verossimilhança
% cai mais retorna rapidamente para valores normais. Isso se deve ao
% alteração do valor de Q que é ajustado para um valor mais adequado. Esta
% técnica é uma versão simples das diferentes maneiras propóstas de
% identificar uma manobra e com isso adaptar parâmetros do modelo e/ou
% filtro de kalman. Aqui pudemos demonstrar alguns benefícios desta
% técnica.

%%
close all