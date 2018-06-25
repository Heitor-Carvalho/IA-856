%% Estimador Yule - Walker

%% Sistema de segunda ordem
% yk + 0.6yk + 0.3yk = wk

% Simulando o sistema
N = 1e4;
wk = randn(N, 1);
yk = filter([1],[1 0.6 0.3], wk);


%% Estimador de AR(1)
Y = convmtx([0; yk], 1);
Y = Y(1:end-1, :);

% Coeficientes AR
theta = Y\yk

% Erro de predição
error = yk - Y*theta;

[autoc, lag] = xcorr(error, error, 10, 'unbiased');

% Verificando se o erro é um ruído branco
subplot(1, 2, 1)
hist(error, 20)
title('Histograma do erro')
grid

subplot(1, 2, 2)
stem(lag, autoc)
title('Autocorrelação')
grid

var(error)

% Aqui tentamos aproximar um sistema AR de segunda orderm por um de primera
% ordem, podemos ver que para neste caso o resíduo não foi perfeitamente branco. 
% e que ainda há correalção entre as amostras do resíduo.

%% Estimador de AR(2)
Y = convmtx([0; yk], 2);
Y = Y(1:end-2, :);

% Coeficientes AR
theta = Y\yk

% Erro de predição
error = yk - Y*theta;

[autoc, lag] = xcorr(error, error, 10, 'unbiased');

% Verificando se o erro é um ruído branco
subplot(1, 2, 1)
hist(error, 20)
title('Histograma do erro')
grid

subplot(1, 2, 2)
stem(lag, autoc)
title('Autocorrelação')
grid

var(error)

% Como pode ser visto no histograma, este se aproxima de uma distribuição
% gaussian, o que indica um ruído branco. Além disso, vemos que
% a variância do erro é 1 e sua função de autocorrelação se aproxima de 
% um impulso para o atraso 0.
