function neg_log_likelihood = motor_biasallp_vb(x, data1, data2)
%MOTOR_BIASALLP_VB Calculates cost for a combined TR+TG model.
%
%   Inputs:
%       x       - A 1x7 vector of model parameters:
%                 x(1)-x(5): Parameters for the TR component (see motor_biasallp).
%                    x(1),x(2): refenrence point
%                    x(3),x(4): Base bias vector.
%                    x(5):     1
%                 x(6)-x(7): Parameters for the TG component.
%                    x(6):     Width/periodicity of the TG function.
%                    x(7):     Amplitude of the TG function.
%       data1   - (Unused) online data
%       data2   - in person data
%
%   Output:
%       neg_log_likelihood 

%% 1. Unpack Parameters
% TR (Transfer) component parameters
v         = [x(1); x(2)];
bias_vec  = [x(3); x(4)];
power_exp = x(5);

% TG (Target Generalization) component parameters
tg_width = x(6);
tg_amp   = x(7);

%% 2. Generate Target Bias

% Define the shape of the bias over 90 degrees.
p1 = 0:tg_width;
p2 = (tg_width-0.01) : (-2*tg_width/(91-2*tg_width-1)) : (-tg_width);
p3 = (-tg_width) : (-0.1); 
single_wave = [p1, p2, p3] * tg_amp;

% Repeat the bias four times to cover the full 360-degree circle.
visual_bias_wave = repmat(single_wave, 1, 4);

% Sample the wave at the 24 target locations to get the specific bias for each.
tg_bias_per_target = visual_bias_wave(1:15:360); % 'vb' in original code

%% 3. Generate Transfer (TR) Bias based on Perceived Targets
% Define the actual angles of the 24 targets.
Ttheta_actual = 0:pi/12:(2*pi - 0.001);

% actual angle plus the TG bias from the previous step.
Ttheta_perceived = Ttheta_actual + tg_bias_per_target;

% Convert the *perceived* target angles to Cartesian coordinates. The TR
% model operates on these biased target representations.
[tx, ty] = pol2cart(Ttheta_perceived, 1);
T_perceived = [tx; ty];

% The rest of the TR model is identical to motor_biasallp, but uses
% the perceived target locations (T_perceived).
start_point_bias = sqrt(nansum(v.^2)).^power_exp .* bias_vec;
distance_to_v = sqrt(nansum((T_perceived - v).^2)).^power_exp;
target_bias = bias_vec .* distance_to_v;
planned_endpoint = T_perceived + target_bias;
planned_movement_vector = planned_endpoint - start_point_bias;

%% 4. Calculate Final Predicted Error
% Convert the final planned movement vector back to an angle.
[Mtheta, ~] = cart2pol(planned_movement_vector(1,:), planned_movement_vector(2,:));

% The final predicted error is the difference between the planned motor
% output angle (Mtheta) and the *actual* target angle (Ttheta_actual).
predicted_error_theta = Mtheta - Ttheta_actual;

% Wrap angles to the range [-pi, pi] 
predicted_error_theta(predicted_error_theta < -pi) = predicted_error_theta(predicted_error_theta < -pi) + 2*pi;
predicted_error_theta(predicted_error_theta > pi)  = predicted_error_theta(predicted_error_theta > pi)  - 2*pi;

%% 5. Log-Likelihood
% Calculate the residuals (difference between prediction and data).
residuals = predicted_error_theta - data2;

% Estimate the noise (sigma) from the residuals.
sigma = nanstd(residuals(:));

% Calculate the negative log-likelihood assuming Gaussian noise.
log_likelihood = log(normpdf(residuals, 0, sigma));
neg_log_likelihood = -nansum(log_likelihood(:));

end