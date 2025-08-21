function neg_log_likelihood = motor_biasallp(x, data1, data2)
%MOTOR_BIASALLP Calculates cost for a TR (Transfer-Only) model.
%
%   Inputs:
%       x       - A 1x5 vector of model parameters:
%                 x(1), x(2): reference point.
%                 x(3), x(4): base bias vector.
%                 x(5):     1 constant
%       data1   - (Unused) online data
%       data2   - In-person data
%
%   Output:
%       neg_log_likelihood

%% 1. Unpack Parameters and Define Targets
% The parameter vector 'x' is unpacked into meaningful variables.
v         = [x(1); x(2)]; % reference point
bias_vec  = [x(3); x(4)]; % Base bias vector
power_exp = x(5);         % Exponent for distance scaling = 1

% Define the 24 standard target locations on a unit circle.
Ttheta = 0:pi/12:(2*pi - 0.001);
[tx, ty] = pol2cart(Ttheta, 1);
T = [tx; ty]; % Target locations in Cartesian coordinates [2 x 24]

%% 2. Generate Model Predictions
% The model assumes bias is applied to the representation of both the
% start and end points of the movement.

% Calculate the bias for the starting point (origin) based on its distance to v.

start_point_bias = sqrt(nansum(v.^2)).^power_exp .* bias_vec;

% Calculate the bias for each target based on its distance to v.

distance_to_v = sqrt(nansum((T - v).^2)).^power_exp;
target_bias = bias_vec .* distance_to_v;

% The planned endpoint is the target location plus the target-specific bias.
planned_endpoint = T + target_bias;

% The final planned movement vector is the difference between the biased
% endpoint and the biased starting point.
planned_movement_vector = planned_endpoint - start_point_bias;

% Convert the planned movement vector back to an angle to find the predicted error.
[Mtheta, ~] = cart2pol(planned_movement_vector(1,:), planned_movement_vector(2,:));
predicted_error_theta = Mtheta - Ttheta;

% Wrap angles to the range [-pi, pi].
predicted_error_theta(predicted_error_theta < -pi) = predicted_error_theta(predicted_error_theta < -pi) + 2*pi;
predicted_error_theta(predicted_error_theta > pi)  = predicted_error_theta(predicted_error_theta > pi)  - 2*pi;

%% 3. Calculate Log-Likelihood
% This step evaluates how well the model's predictions match the observed data.

% Calculate the residuals (difference between prediction and data).
residuals = predicted_error_theta - data2;

% Estimate the noise.
sigma = nanstd(residuals(:));

% Calculate the log-likelihood
log_likelihood = log(normpdf(residuals, 0, sigma));

% The final cost is the sum of the negative log-likelihoods.
neg_log_likelihood = -nansum(log_likelihood(:));

end


