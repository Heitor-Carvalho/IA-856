%% Generating simple seismic trace
wav = mexihat(-30, 30, 100);           % Wave format

% To do: Generate some distorte wave format to test

imp = [1 0 0 2 0 0 2 0 0 2];           % Perfect impulse response
impp = upsample(imp, 50);              % Upsampled impulse response
trace = conv(impp, wav);               % Perfect trace

imp_att = [1 0 0 0.8 0 0 0.1 0 0 0.05]; % Attenuated impulse response
impp = upsample(imp_att, 50);          % Upsampled impulse response
trace_att = conv(impp, wav);           % Attenuated trace

imp = [1 0 0 1 1 0 1 1 0 1];           % Perfect impulse response
impp = upsample(imp, 50);              % Upsampled impulse response
trace_one_prim = conv(impp, wav);      % Perfect trace

imp = [1 0 0 0.8 1 0 0.4 0.6 0 0.1];   % Perfect impulse response
impp = upsample(imp, 50);              % Upsampled impulse response
trace_one_prim_att = conv(impp, wav);  % Perfect trace


filter_len = 1;
prediction_step = 75;

%% Perfect trace case
% Organizing data for prediction
M = convmtx(trace', filter_len);
MM_perf = M(1:end-prediction_step-filter_len+1, :);

% Calculating linear regression coeficients
gain_perf = inv(MM_perf'*MM_perf)*MM_perf'*trace(1+prediction_step:end)';
pred_error_perf = MM_perf*gain_perf-trace(1+prediction_step:end)';
pred_error_perf = trace' - [zeros(prediction_step, 1); MM_perf*gain_perf];
pred_error_perf = pred_error_perf(1:end-100);

figure(1)
plot(pred_error_perf)
ylim([-1 1])
grid

if filter_len == 2
  figure(111)
  xy_range = [min(min(MM_perf(:, 1), MM_perf(:, 2))), max(max(MM_perf(:, 1), MM_perf(:, 2)))];
  [xx, yy] = meshgrid(xy_range(1):(xy_range(2)-xy_range(1))/100:xy_range(2));
  filter_surf = [xx(:) yy(:)]*gain_perf;
  filter_surf = reshape(filter_surf, size(xx));
  mesh(xx, yy, filter_surf)
  hold on
  plot3(MM_perf(:, 1), MM_perf(:, 2), trace(1+prediction_step:end)')
  grid
end

%% Attenuated trace case
% Organizing data for prediction
M = convmtx(trace_att', filter_len);
MM_att = M(1:end-prediction_step-filter_len+1, :);

% Calculating linear regression coeficients
gain_att = inv(MM_att'*MM_att)*MM_att'*trace_att(1+prediction_step:end)';
pred_error_att = trace_att' - [zeros(prediction_step, 1); MM*gain_att];
pred_error_att = pred_error_att(1:end-100);

figure(2)
plot(pred_error_att)
ylim([-1 1])
grid

if filter_len == 2
  figure(112)
  xy_range = [min(min(MM_att(:, 1), MM_att(:, 2))), max(max(MM_att(:, 1), MM_att(:, 2)))];
  [xx, yy] = meshgrid(xy_range(1):(xy_range(2)-xy_range(1))/100:xy_range(2));
  filter_surf = [xx(:) yy(:)]*gain_perf;
  filter_surf = reshape(filter_surf, size(xx));
  mesh(xx, yy, filter_surf)
  hold on
  plot3(MM_att(:, 1), MM_att(:, 2), trace_att(1+prediction_step:end)')
  grid
end

%% One primary trace case
% Organizing data for prediction
M = convmtx(trace_one_prim', filter_len);
MM_one_prim = M(1:end-prediction_step-filter_len+1, :);

% Calculating linear regression coeficients
gain_prim = inv(MM_one_prim'*MM_one_prim)*MM_one_prim'*trace_one_prim(1+prediction_step:end)';
pred_error_prim = trace_one_prim' - [zeros(prediction_step, 1); MM_one_prim*gain_prim];
pred_error_prim = pred_error_prim(1:end-100);

figure(3)
plot(pred_error_prim)
ylim([-1 1])
grid

if filter_len == 2
  figure(113)
  xy_range = [min(min(MM_one_prim(:, 1), MM_one_prim(:, 2))), max(max(MM_one_prim(:, 1), MM_one_prim(:, 2)))];
  [xx, yy] = meshgrid(xy_range(1):(xy_range(2)-xy_range(1))/100:xy_range(2));
  filter_surf = [xx(:) yy(:)]*gain_prim;
  filter_surf = reshape(filter_surf, size(xx));
  mesh(xx, yy, filter_surf)
  hold on
  plot3(MM_one_prim(:, 1), MM_one_prim(:, 2), trace_one_prim(1+prediction_step:end)')
  grid
end


%% One primary trace att case
% Organizing data for prediction
M = convmtx(trace_one_prim_att', filter_len);
MM_one_prim_att = M(1:end-prediction_step-filter_len+1, :);

% Calculating linear regression coeficients
gain_prim_att = inv(MM_one_prim_att'*MM_one_prim_att)*MM_one_prim_att'*trace_one_prim_att(1+prediction_step:end)';
pred_error_prim_att = trace_one_prim_att' - [zeros(prediction_step, 1); MM_one_prim_att*gain_prim_att];
pred_error_prim_att = pred_error_prim_att(1:end-100);

figure(3)
plot(pred_error_prim_att)
ylim([-1 1])
grid

if filter_len == 2
  figure(113)
  xy_range = [min(min(MM_one_prim(:, 1), MM_one_prim(:, 2))), max(max(MM_one_prim(:, 1), MM_one_prim(:, 2)))];
  [xx, yy] = meshgrid(xy_range(1):(xy_range(2)-xy_range(1))/100:xy_range(2));
  filter_surf = [xx(:) yy(:)]*gain_prim;
  filter_surf = reshape(filter_surf, size(xx));
  mesh(xx, yy, filter_surf)
  hold on
  plot3(MM_one_prim(:, 1), MM_one_prim(:, 2), trace_one_prim(1+prediction_step:end)')
  grid
end

%% Conclustion: even in the perfect scenario a linear regression cannot catch up
% to the attenuation. But maybe a almost linear regression can!
% When we do a linear regression we are assuming we are obserarving a static system,
% the idea here is to consider that we are observing a very slow dynamical system
% this way we will be able to better adapt to the attenuation process.
% Our dynamical system must be very slow, since we do not want to catch up with primary
% present in our data, this primary represent a very suddenly change, that is why a linear
% regression do not kill this part of the signal. In this first framework our primary will
% be fully contained in our prediction error.

% The idea in the future is to desgin another model able to follow this non predictable
% primary and let the Interacting Multiple Model to detect which model is correct
% This fisrt model, is a very low dynamical system, we are triying to get a very
% close to linear regression result. The second model is free to follow some changes


%% Filtering perfect case with Kalman, RLS like
x0 = zeros(filter_len, 1);
P0 = eye(length(x0))*1e2;
F = eye(length(x0));
H = eye(length(x0));
Q = eye(length(x0))*0;                 % We are assuming a perferc model
R = eye(length(x0))*1e-3;              % Very low noisy measurments

Hobs = zeros(length(x0), length(x0), size(MM_perf, 1));
for i = 1:size(Hobs, 3)
  Hobs(1,:,i) = MM_perf(i, :);
end

[x_est, likelyhood, residual_hist, ~, P] = kalman_filter_wrapper_dyn_h(x0, P0, trace(1+prediction_step:end)', F, Hobs, Q, R);

y = zeros(1, size(x_est, 2));
for i = 1:size(Hobs, 3)
  y(i) = Hobs(1,:,i)*x_est(:, i);
end

pred_error_perf_kalman = trace' - [zeros(prediction_step, 1); y'];



Hobs = zeros(length(x0), length(x0), size(MM_att, 1));
for i = 1:size(Hobs, 3)
  Hobs(1,:,i) = MM_att(i, :);
end

[x_est, likelyhood, residual_hist, ~, P] = kalman_filter_wrapper_dyn_h(x0, P0, trace_att(1+prediction_step:end)', F, Hobs, Q, R);

y = zeros(1, size(x_est, 2));
for i = 1:size(Hobs, 3)
  y(i) = Hobs(1,:,i)*x_est(:, i);
end

pred_error_att_kalman = trace_att' - [zeros(prediction_step, 1); y'];




Hobs = zeros(length(x0), length(x0), size(MM_one_prim, 1));
for i = 1:size(Hobs, 3)
  Hobs(1,:,i) = MM_one_prim(i, :);
end

[x_est, likelyhood, residual_hist, ~, P] = kalman_filter_wrapper_dyn_h(x0, P0, trace_one_prim(1+prediction_step:end)', F, Hobs, Q, R);

y = zeros(1, size(x_est, 2));
for i = 1:size(Hobs, 3)
  y(i) = Hobs(1,:,i)*x_est(:, i);
end

pred_error_one_prim_kalman = trace_one_prim' - [zeros(prediction_step, 1); y'];




Hobs = zeros(length(x0), length(x0), size(MM_one_prim_att, 1));
for i = 1:size(Hobs, 3)
  Hobs(1,:,i) = MM_one_prim_att(i, :);
end

[x_est, likelyhood, residual_hist, ~, P] = kalman_filter_wrapper_dyn_h(x0, P0, trace_one_prim_att(1+prediction_step:end)', F, Hobs, Q, R);

y = zeros(1, size(x_est, 2));
for i = 1:size(Hobs, 3)
  y(i) = Hobs(1,:,i)*x_est(:, i);
end

pred_error_one_prim_att_kalman = trace_one_prim_att' - [zeros(prediction_step, 1); y'];
