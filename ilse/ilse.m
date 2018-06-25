%% Execício 3-15 
clear all
clc

%% Parâmetros cenário 1:
xy1 = [20000, 20000];                         % Posição real
xy_sc1_mp = [0, 0; 15000, 5000; 30000, 0];    % Posição dos pontos de medidas

%% Calculando Jacobian/Gradient matrix
syms x x0 y y0

r(x,y,x0,y0) = sqrt((x-x0)^2 + (y-y0)^2)

drdx = diff(r,x);
drdy = diff(r,y);

J = [drdx; drdy]

%% Encontrando as derivadas usadnas no cálcula da matriz Hessiana
drdxx = diff(drdx, x);
drdyx = diff(drdy, x);
drdxy = diff(drdx, y);
drdyy = diff(drdy, y);

%% Carregando os parâmetros para o cenario 1
xy_mp = xy_sc1_mp;
xy_real = repmat(xy1, length(xy_mp), 1);

%% Obtendo medida ruidosa
sigma = 40^2;
z = double(r(xy_real(:,1), xy_real(:,2), xy_mp(:,1), xy_mp(:,2))) + sqrt(sigma)*randn(3, 1);

%% Caso 1 - Iterative least square
N = 10;
x0 = [20100 19910];
x_est = x0;

R = sigma^2*eye(3);

x_est_hist = zeros(N, length(x0));
for i=1:N
  x_est_hist(i,:) = x_est;
  drdx_num = double(drdx(x_est(:,1), x_est(:,2), xy_mp(:,1), xy_mp(:,2)));
  drdy_num = double(drdy(x_est(:,1), x_est(:,2), xy_mp(:,1), xy_mp(:,2)));
  JJ = [drdx_num, drdy_num];
  x_est_stack = repmat(x_est, length(xy_mp), 1);
  r_est = double(r(x_est(:,1), x_est(:,2), xy_mp(:,1), xy_mp(:,2)));
  x_est = x_est + (inv(JJ'*inv(R)*JJ)*JJ'*inv(R)*(z-r_est))';
end

figure(1)
plot(x_est_hist, '--o')
legend('Estados estimados')
grid

%% Erro do valor estimado IRLS
x_est - xy_real(1, :)

%% Caso 2 - Estimador ML with Newton-Raphson
N = 10;
x0 = [20100 19910];
x_est = x0;
sigma = 40^2;

x_est_hist = zeros(N, length(x0));
for i = 1:N
  x_est_hist(i,:) = x_est;
  r_est = double(r(x_est(:,1), x_est(:,2), xy_mp(:,1), xy_mp(:,2)));
  error = z - r_est;

  dldxx_num = error.*double(drdxx(x_est(:,1), x_est(:,2), xy_mp(:,1), xy_mp(:,2))) - double(drdx(x_est(:,1), x_est(:,2), xy_mp(:,1), xy_mp(:,2))).^2;
  dldxx_num = sum(dldxx_num)/sigma.^2;
  dldyy_num = error.*double(drdyy(x_est(:,1), x_est(:,2), xy_mp(:,1), xy_mp(:,1))) - double(drdy(x_est(:,1), x_est(:,2), xy_mp(:,1), xy_mp(:,2))).^2;
  dldyy_num = sum(dldyy_num)/sigma.^2;
  dldxy_num = error.*double(drdxy(x_est(:,1), x_est(:,2), xy_mp(:,1), xy_mp(:,2))) -  ...
              double(drdx(x_est(:,1), x_est(:,2), xy_mp(:,1), xy_mp(:,1))).*double(drdy(x_est(:,1), x_est(:,2), xy_mp(:,1), xy_mp(:,2)));
  dldxy_num = sum(dldxy_num)/sigma.^2;
  dldyx_num = dldxy_num;

  dldx_num = error.*double(drdx(x_est(:,1), x_est(:,2), xy_mp(:,1), xy_mp(:,2)));
  dldx_num = sum(dldx_num)/sigma.^2;
  dldy_num = error.*double(drdy(x_est(:,1), x_est(:,2), xy_mp(:,1), xy_mp(:,2)));
  dldy_num = sum(dldy_num)/sigma.^2;

  hessian = [dldxx_num dldxy_num; ...
             dldyx_num dldyy_num];
  gradiant = [dldx_num; dldy_num];

  x_est = x_est - (inv(hessian)*gradiant)';
end

figure(2)
plot(x_est_hist, '--o')
legend('Estados estimados')
grid

%% Erro do valor estimado ML
x_est - xy_real(1, :)

% Neste cenário, podemos ver que ambos os métodos convergiram rapidamente e
% tiveram comportamentos semelhantes. Além disso, os erros de estimação
% foram os mesmos.

%% Calculando o Cramer-Rao Lower Bound
drdx_num = double(drdx(xy_real(:,1), xy_real(:,2), xy_mp(:,1), xy_mp(:,2)));
drdy_num = double(drdy(xy_real(:,1), xy_real(:,2), xy_mp(:,1), xy_mp(:,2)));
JJ = [drdx_num, drdy_num];

CRLB = zeros(size(xy_real, 2));
for i = 1:length(xy_mp)
  CRLB = CRLB + JJ(i,:)'*JJ(i,:)/sigma^2; 
end
CRLB = inv(CRLB)

% Calculando Sigma CRLB
sigma_crlb_x = sqrt(CRLB(1,1))
sigma_crlb_y = sqrt(CRLB(2,2))
