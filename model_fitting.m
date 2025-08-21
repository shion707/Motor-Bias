%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This script fits a suite of seven different models to the group-averaged
%  in-person dataset to explain motor biases.
%  Each model was fit for 10 times in this demo
%  A robust fit reqires 150-200 fits!
%  The script is organized as follows:
%  1. Load and preprocess all required datasets (In-Person and Online).
%  2. Create a summary figure showing the preprocessed data for each set.
%  3. Sequentially fit each of the 7 models to the in-person data.
%  4. Plot the results for each model fit in a separate figure.
%  5. Calculate the Bayesian Information Criterion (BIC) for each model.

clear; close all; clc;
%% =========================================================================
%  1. Load and Preprocess All Datasets
%  =========================================================================
% --- In-Person 8-target data ---
T = readtable('inperson_8target.csv');
trial_num = 80; subnum_ip8 = 150;
hand_common_angle = T.Hand/180*pi;
hand_common_angle = reshape(hand_common_angle, trial_num, subnum_ip8);
hand_common_angle = hand_common_angle(1:trial_num/2, :);
T_angle = T.ti; T_angle = reshape(T_angle, trial_num, subnum_ip8);
T_angle = T_angle(1:trial_num/2, :); T_angle = T_angle(:);
[~, j] = sort(T_angle);
hand_common_angle = hand_common_angle(j);
hand_common_angle = reshape(hand_common_angle, 12000/8/2, 8);
mhand = nanmean(hand_common_angle); shand = nanstd(hand_common_angle);
hand_common_angle(hand_common_angle > mhand + 3*shand) = nan;
hand_common_angle(hand_common_angle < mhand - 3*shand) = nan;
hand_common_angle = reshape(hand_common_angle, 12000/8/2/subnum_ip8, subnum_ip8, 8);
hand_common_angle_p8 = squeeze(nanmean(hand_common_angle));   % [Nsubj x 8]
% --- In-Person 24-target data ---
T = readtable('inperson_24target.csv');
trial_num = 96 * 2; subnum_ip24 = 56;
hand_common_angle = T.Hand/180*pi;
hand_common_angle = reshape(hand_common_angle, trial_num, subnum_ip24);
hand_common_angle = hand_common_angle(1:trial_num/2, :);
T_angle = T.ti; T_angle = reshape(T_angle, trial_num, subnum_ip24);
T_angle = T_angle(1:trial_num/2, :); T_angle = T_angle(:);
[~, j] = sort(T_angle);
hand_common_angle = hand_common_angle(j);
hand_common_angle = reshape(hand_common_angle, 10752/24/2, 24);
mhand = nanmean(hand_common_angle); shand = nanstd(hand_common_angle);
hand_common_angle(hand_common_angle > mhand + 3*shand) = nan;
hand_common_angle(hand_common_angle < mhand - 3*shand) = nan;
hand_common_angle = reshape(hand_common_angle, 10752/24/2/subnum_ip24, subnum_ip24, 24);
hand_common_angle_p24 = squeeze(nanmean(hand_common_angle));  % [Nsubj x 24]
mhand_p24 = nanmedian(hand_common_angle_p24);
% --- Online 8-target data (loaded but not used for fitting) ---
T = readtable('online_8target.csv');
trial_num = 320; subnum_on8 = 221;
hand_common_angle = T.Hand/180*pi;
hand_common_angle = reshape(hand_common_angle, trial_num, subnum_on8);
hand_common_angle = hand_common_angle(1:trial_num/2, :);
T_angle = T.ti;
T_angle = reshape(T_angle, trial_num, subnum_on8);
T_angle = T_angle(1:trial_num/2, :); T_angle = T_angle(:);
[~, j] = sort(T_angle);
hand_common_angle = hand_common_angle(j);
hand_common_angle = reshape(hand_common_angle, 70720/8/2, 8);
mhand = nanmean(hand_common_angle); shand = nanstd(hand_common_angle);
hand_common_angle(hand_common_angle > mhand + 3*shand) = nan;
hand_common_angle(hand_common_angle < mhand - 3*shand) = nan;
hand_common_angle = reshape(hand_common_angle, 70720/8/2/subnum_on8, subnum_on8, 8);
hand_common_angle_o8 = squeeze(nanmean(hand_common_angle));   % [Nsubj x 8]
% --- Online 24-target data (loaded but not used for fitting) ---
T = readtable('online_24target.csv');
trial_num = 960; subnum_on24 = 69;
hand_common_angle = T.Hand/180*pi;
hand_common_angle = reshape(hand_common_angle, trial_num, subnum_on24);
hand_common_angle = hand_common_angle(1:trial_num/2, :);
T_angle = T.ti;
T_angle = reshape(T_angle, trial_num, subnum_on24);
T_angle = T_angle(1:trial_num/2, :); T_angle = T_angle(:);
[~, j] = sort(T_angle);
hand_common_angle = hand_common_angle(j);
hand_common_angle = reshape(hand_common_angle, 66240/24/2/subnum_on24, subnum_on24, 24);
hand_common_angle_o24 = squeeze(nanmedian(hand_common_angle)); % [Nsubj x 24]
mhand_o24 = nanmedian(hand_common_angle_o24);
%% =========================================================================
%  1b. Visualize All Datasets Before Fitting
%  =========================================================================
Ttheta8  = 0:pi/4:2*pi-0.001;  x8_deg  = Ttheta8/pi*180;
Ttheta24 = 0:pi/12:2*pi-0.001; x24_deg = Ttheta24/pi*180;
% In-Person 8-target
m_ip8  = nanmedian(hand_common_angle_p8);  % 1x8
se_ip8 = nanstd(hand_common_angle_p8) ./ sqrt(size(hand_common_angle_p8,1));
figure(1); clf; hold on; title('In-Person 8-target (pre-fit)');
errorbar(x8_deg, m_ip8/pi*180, se_ip8/pi*180, 'o', ...
    'MarkerFaceColor','m', 'Color','m', 'CapSize',0, 'LineStyle','none', 'LineWidth', 1);
yline(0,'-k','LineWidth',0.5); xlim([-45 360]); ylim([-20 20]);
xlabel('Target Angle (°)'); ylabel('Bias (°)'); grid on; box on;
legend('Median \pm SE','Location','northwest');
text(5, 18, sprintf('n = %d', size(hand_common_angle_p8,1)));
% In-Person 24-target
m_ip24  = nanmedian(hand_common_angle_p24);  % 1x24
se_ip24 = nanstd(hand_common_angle_p24) ./ sqrt(size(hand_common_angle_p24,1));
figure(2); clf; hold on; title('In-Person 24-target (pre-fit)');
errorbar(x24_deg, m_ip24/pi*180, se_ip24/pi*180, 'o', ...
    'MarkerFaceColor','m', 'Color','m', 'CapSize',0, 'LineStyle','none', 'LineWidth', 1);
yline(0,'-k','LineWidth',0.5); xlim([-45 360]); ylim([-20 20]);
xlabel('Target Angle (°)'); ylabel('Bias (°)'); grid on; box on;
legend('Median \pm SE','Location','northwest');
text(5, 18, sprintf('n = %d', size(hand_common_angle_p24,1)));
% Online 8-target
m_on8  = nanmedian(hand_common_angle_o8);
se_on8 = nanstd(hand_common_angle_o8) ./ sqrt(size(hand_common_angle_o8,1));
figure(3); clf; hold on; title('Online 8-target (pre-fit)');
errorbar(x8_deg, m_on8/pi*180, se_on8/pi*180, 'o', ...
    'MarkerFaceColor','m', 'Color','m', 'CapSize',0, 'LineStyle','none', 'LineWidth', 1);
yline(0,'-k','LineWidth',0.5); xlim([-45 360]); ylim([-20 20]);
xlabel('Target Angle (°)'); ylabel('Bias (°)'); grid on; box on;
legend('Median \pm SE','Location','northwest');
text(5, 18, sprintf('n = %d', size(hand_common_angle_o8,1)));
% Online 24-target
m_on24  = nanmedian(hand_common_angle_o24);
se_on24 = nanstd(hand_common_angle_o24) ./ sqrt(size(hand_common_angle_o24,1));
figure(4); clf; hold on; title('Online 24-target (pre-fit)');
errorbar(x24_deg, m_on24/pi*180, se_on24/pi*180, 'o', ...
    'MarkerFaceColor','m', 'Color','m', 'CapSize',0, 'LineStyle','none', 'LineWidth', 1);
yline(0,'-k','LineWidth',0.5); xlim([-45 360]); ylim([-20 20]);
xlabel('Target Angle (°)'); ylabel('Bias (°)'); grid on; box on;
legend('Median \pm SE','Location','northwest');
text(5, 18, sprintf('n = %d', size(hand_common_angle_o24,1)));
%% =========================================================================
%  2. Initialize Common Variables (for model fitting)
%  =========================================================================
options = optimset('MaxFunEvals', 10000, 'MaxIter', 50000);
Ttheta24 = 0:pi/12:2*pi-0.001;
n_data_points = length(hand_common_angle_p24(:));
data_to_fit = hand_common_angle_p24;                  % [Nsubj x 24]
Nsubj_fit   = size(data_to_fit,1);
x24_deg     = Ttheta24/pi*180;
mdata_deg   = mhand_p24/pi*180;
se_fit_deg  = (nanstd(data_to_fit) ./ sqrt(Nsubj_fit)) / pi * 180;
%% Helper to add a BIC annotation (top-left in axis coords)
addBIC = @(bicVal) text(min(xlim)+5, max(ylim)-2, sprintf('BIC = %.1f', bicVal), ...
                        'HorizontalAlignment','left','VerticalAlignment','top', ...
                        'BackgroundColor','w','Margin',4);
%% =========================================================================
%  3. Fit & Plot TR Model
%  =========================================================================
disp('Fitting Model 1: TR (Transfer-Only)...');
min_neg_log_likelihood = 1e10;
for i = 1:10
    fun = @(x) motor_biasallp(x, hand_common_angle_p8, data_to_fit);
    x0 = [8, -0.5, -0.0001, -0.00005, 1];
    x1 = [0, -10, -inf, -inf, 1];
    x2 = [15, 0, 0, 0, 1];
    [x, fval] = fminsearchbnd(fun, x0, x1, x2, options);
    if fval < min_neg_log_likelihood
        min_neg_log_likelihood = fval; best_fit_params_tr = x;
    end
end
BIC_TR = 2 * min_neg_log_likelihood * 10 + 4 * log(n_data_points * 10);
x = best_fit_params_tr; v = [x(1); x(2)];
[tx24, ty24] = pol2cart(Ttheta24, 1); T24 = [tx24; ty24];
bias_vec = [x(3); x(4)];
distance24 = sqrt(nansum((T24 - v).^2)).^x(5);
bias24 = bias_vec .* distance24; tp24 = T24 + bias24;
startp = sqrt(nansum(v.^2)).^x(5) .* bias_vec;
planp24 = tp24 - startp;
[Mtheta24, ~] = cart2pol(planp24(1,:), planp24(2,:));
error_theta = Mtheta24 - Ttheta24;
error_theta(error_theta < -pi) = error_theta(error_theta < -pi) + 2*pi;
figure(101); clf; hold on; title('Model 1: TR (Transfer-Only)');
errorbar(x24_deg, mdata_deg, se_fit_deg, 'o', ...
    'MarkerFaceColor',[1,0.5,1], 'Color',[1,0.5,1], 'CapSize',0, 'LineStyle','none', 'LineWidth', 1);
plot(x24_deg, error_theta/pi*180, '-', 'linewidth', 2, 'color', [0.5,0.5,1]);
plot([-45 360], [0 0], '-k', 'LineWidth', 0.5);
axis([-45 360 -15 25]); ylabel('Bias (°)'); xlabel('Target Angle (°)');
set(gca, 'xtick', 0:90:360, 'ytick', -10:10:20, 'LineWidth', 1); box on; grid on;
legend('Median \pm SE','Model fit','Location','northwest');
addBIC(BIC_TR);
%% =========================================================================
%  4. Fit & Plot TR+TG Model
%  =========================================================================
disp('Fitting Model 2: TR+TG (Transfer + Target Generalization)...');
min_neg_log_likelihood = 1e10;
for i = 1:10
    fun = @(x) motor_biasallp_vb(x, hand_common_angle_p8, data_to_fit);
    x0 = [randn*10, randn*2, -randn, -randn, 1, randi(10)+4, rand];
    x1 = [0, -10, -inf, -inf, 1, 4, 0];
    x2 = [15, 0, 0, 0, 1, inf, inf];
    [x, fval] = fminsearchbnd(fun, x0, x1, x2, options);
    if fval < min_neg_log_likelihood
        min_neg_log_likelihood = fval; best_fit_params_tr_tg = x;
    end
end
BIC_TR_TG = 2 * min_neg_log_likelihood * 10 + 6 * log(n_data_points * 10);
x = best_fit_params_tr_tg; v = [x(1); x(2)]; bias_vec = [x(3); x(4)];
p1 = 0:x(6); p2 = (x(6)-0.01):(-2*x(6)/(90-2*x(6)-1)):(-x(6)); p3 = (-x(6)):(-0.1);
target_gen_bias = [p1, p2, p3] * x(7); target_gen_bias = repmat(target_gen_bias, 1, 4);
tg_24_targets = target_gen_bias(1:15:360);
distance24 = sqrt(nansum((T24 - v).^2)).^x(5);
bias24 = bias_vec .* distance24; tp24 = T24 + bias24;
startp = sqrt(nansum(v.^2)).^x(5) .* bias_vec;
planp24 = tp24 - startp;
[Mtheta24, ~] = cart2pol(planp24(1,:), planp24(2,:));
error_theta = Mtheta24 - Ttheta24 + tg_24_targets;
error_theta(error_theta < -pi) = error_theta(error_theta < -pi) + 2*pi;
figure(102); clf; hold on; title('Model 2: TR+TG');
errorbar(x24_deg, mdata_deg, se_fit_deg, 'o', ...
    'MarkerFaceColor',[1,0.5,1], 'Color',[1,0.5,1], 'CapSize',0, 'LineStyle','none', 'LineWidth', 1);
plot(x24_deg, error_theta/pi*180, '-', 'linewidth', 2, 'color', [0.5,0.5,1]);
plot([-45 360], [0 0], '-k', 'LineWidth', 0.5);
axis([-45 360 -15 25]); ylabel('Bias (°)'); xlabel('Target Angle (°)');
set(gca, 'xtick', 0:90:360, 'ytick', -10:10:20, 'LineWidth', 1); box on; grid on;
legend('Median \pm SE','Model fit','Location','northwest');
addBIC(BIC_TR_TG);
%% =========================================================================
%  5. Fit & Plot TG-Only Model
%  =========================================================================
disp('Fitting Model 3: TG-Only (Target Generalization Only)...');
min_neg_log_likelihood = 1e10;
for i = 1:10
    fun = @(x) motor_biasallp_vb_only(x, data_to_fit);
    x0 = [randi(10)+4, rand]; x1 = [4, 0]; x2 = [inf, inf];
    [x, fval] = fminsearchbnd(fun, x0, x1, x2, options);
    if fval < min_neg_log_likelihood
        min_neg_log_likelihood = fval; best_fit_params_tg = x;
    end
end
BIC_TG = 2 * min_neg_log_likelihood * 10 + 2 * log(n_data_points * 10);
x = best_fit_params_tg;
p1 = 0:x(1); p2 = (x(1)-0.01):(-2*x(1)/(91-2*x(1)-1)):(-x(1)); p3 = (-x(1)):(-0.1);
target_gen_bias = [p1, p2, p3] * x(2);
target_gen_bias = repmat(target_gen_bias, 1, 4);
error_theta = target_gen_bias(1:15:360);
figure(103); clf; hold on; title('Model 3: TG-Only');
errorbar(x24_deg, mdata_deg, se_fit_deg, 'o', ...
    'MarkerFaceColor',[1,0.5,1], 'Color',[1,0.5,1], 'CapSize',0, 'LineStyle','none', 'LineWidth', 1);
plot(x24_deg, error_theta/pi*180, '-', 'linewidth', 2, 'color', [0.5,0.5,1]);
plot([-45 360], [0 0], '-k', 'LineWidth', 0.5);
axis([-45 360 -20 20]); ylabel('Bias (°)'); xlabel('Target Angle (°)');
set(gca, 'xtick', 0:90:360, 'ytick', -20:10:20, 'LineWidth', 1); box on; grid on;
legend('Median \pm SE','Model fit','Location','northwest');
addBIC(BIC_TG);
%% =========================================================================
%  6. Fit & Plot Proprioceptive-Only Model
%  =========================================================================
disp('Fitting Model 4: Proprioceptive-Only...');
min_neg_log_likelihood = 1e10;
for i = 1:100
    fun = @(x) motor_biasprop(data_to_fit, x);
    x0 = [0.1*randn, 0.1*randn 0 0 0 0]; x1 = [-0.9, -0.9 0 0 0 0]; x2 = [0.9, 0.9 0 0 0 0];
    [x, fval] = fminsearchbnd(fun, x0, x1, x2, options);
    if fval < min_neg_log_likelihood
        min_neg_log_likelihood = fval; best_fit_params_prop = x;
    end
end
BIC_Prop = 2 * min_neg_log_likelihood * 10 + 2 * log(n_data_points * 10);
x = best_fit_params_prop;
startp = [x(1); x(2)];
planp = T24 - startp;
[ttemp, ~] = cart2pol(planp(1,:), planp(2,:));
error_theta = ttemp - Ttheta24;
error_theta(error_theta < -pi) = error_theta(error_theta < -pi) + 2*pi;
figure(104); clf; hold on; title('Model 4: Proprioceptive-Only');
errorbar(x24_deg, mdata_deg, se_fit_deg, 'o', ...
    'MarkerFaceColor',[1,0.5,1], 'Color',[1,0.5,1], 'CapSize',0, 'LineStyle','none', 'LineWidth', 1);
plot(x24_deg, error_theta/pi*180, '-', 'linewidth', 2, 'color', [0.5,0.5,1]);
plot([-45 360], [0 0], '-k', 'LineWidth', 0.5);
axis([-45 360 -15 15]); ylabel('Bias (°)'); xlabel('Target Angle (°)');
set(gca, 'xtick', 0:90:360, 'ytick', -15:5:15, 'LineWidth', 1); box on; grid on;
legend('Median \pm SE','Model fit','Location','northwest');
addBIC(BIC_Prop);
%% =========================================================================
%  7. Fit & Plot Proprioceptive+TG Model
%  =========================================================================
disp('Fitting Model 5: Proprioceptive+TG...');
min_neg_log_likelihood = 1e10;
for i = 1:10
    fun = @(x) motor_biasprop_VB(data_to_fit, x);
    x0 = [0.1*randn, 0.1*randn, 0, 0, 4+20*rand, 0.0035];
    x1 = [-0.9, -0.9, 0, 0, 4, 0]; x2 = [0.9, 0.9, 0, 0, inf, inf];
    [x, fval] = fminsearchbnd(fun, x0, x1, x2, options);
    if fval < min_neg_log_likelihood
        min_neg_log_likelihood = fval; best_fit_params_prop_tg = x;
    end
end
BIC_Prop_TG = 2 * min_neg_log_likelihood * 10 + 4 * log(n_data_points * 10);
x = best_fit_params_prop_tg;
startp = [x(1); x(2)]; planp = T24 - startp;
[ttemp, rtemp] = cart2pol(planp(1,:), planp(2,:));
rtemp = x(3)*rtemp; ttemp = ttemp + x(4) + rtemp;
p1 = 0:x(5); p2 = (x(5)-0.01):(-2*x(5)/(91-2*x(5)-1)):(-x(5)); p3 = (-x(5)):(-0.1);
target_gen_bias = [p1, p2, p3] * x(6); target_gen_bias = repmat(target_gen_bias, 1, 4);
tg_24_targets = target_gen_bias(1:15:360);
error_theta = ttemp - Ttheta24 + tg_24_targets;
error_theta(error_theta < -pi) = error_theta(error_theta < -pi) + 2*pi;
figure(105); clf; hold on; title('Model 5: Proprioceptive+TG');
errorbar(x24_deg, mdata_deg, se_fit_deg, 'o', ...
    'MarkerFaceColor',[1,0.5,1], 'Color',[1,0.5,1], 'CapSize',0, 'LineStyle','none', 'LineWidth', 1);
plot(x24_deg, error_theta/pi*180, '-', 'linewidth', 2, 'color', [0.5,0.5,1]);
plot([-45 360], [0 0], '-k', 'LineWidth', 0.5);
axis([-45 360 -15 15]); ylabel('Bias (°)'); xlabel('Target Angle (°)');
set(gca, 'xtick', 0:90:360, 'ytick', -15:5:15, 'LineWidth', 1); box on; grid on;
legend('Median \pm SE','Model fit','Location','northwest');
addBIC(BIC_Prop_TG);
%% =========================================================================
%  8. Fit & Plot Sober Model
%  =========================================================================
disp('Fitting Model 6: Sober (Kinematic)...');
min_neg_log_likelihood = 1e10;
for i = 1:10
    fun = @(x) motor_biasSober2(data_to_fit, x);
    x0 = [pi/2*rand, pi*9/10*rand, 0.1*sign(randn), 0.1*sign(randn)];
    x1 = [0, pi/2, -1, -1]; x2 = [pi/2, pi*9/10, 1, 1];
    [x, fval] = fminsearchbnd(fun, x0, x1, x2, options);
    if fval < min_neg_log_likelihood
        min_neg_log_likelihood = fval; best_fit_params_sober = x;
    end
end
BIC_Sober = 2 * min_neg_log_likelihood * 10 + 4 * log(n_data_points * 10);
x = best_fit_params_sober;
l1 = 3; l2 = 3; a = x(1); b = x(2); biasa = x(3); biasb = x(4);
x_shoulder = cos(a)*l1 + cos(b)*l2; y_shoulder = sin(a)*l1 + sin(b)*l2;
[tx, ty] = pol2cart(Ttheta24, 1);
error_theta = zeros(1, 24);
for ta = 1:24
    x_target = x_shoulder + tx(ta); y_target = y_shoulder + ty(ta);
    d = sqrt(x_target^2 + y_target^2);
    theta1 = acos((l1^2 - l2^2 + d^2) / (2*l1*d));
    theta2 = atan2(y_target, x_target);
    a2 = theta2 - theta1; da = a2 - a;
    theta3 = acos((l1^2 + l2^2 - d^2) / (2*l1*l2));
    b2 = pi - theta3 + a2; db = b2 - b;
    rx = cos(a+da+biasa)*l1 + cos(b+db+biasb)*l2;
    ry = sin(a+da+biasa)*l1 + sin(b+db+biasb)*l2;
    [etheta, ~] = cart2pol(rx - x_shoulder, ry - y_shoulder);
    error_theta(ta) = etheta - Ttheta24(ta);
end
error_theta(error_theta < -pi) = error_theta(error_theta < -pi) + 2*pi;
figure(106); clf; hold on; title('Model 6: Sober (Kinematic)');
errorbar(x24_deg, mdata_deg, se_fit_deg, 'o', ...
    'MarkerFaceColor',[1,0.5,1], 'Color',[1,0.5,1], 'CapSize',0, 'LineStyle','none', 'LineWidth', 1);
plot(x24_deg, error_theta/pi*180, '-', 'linewidth', 2, 'color', [0.5,0.5,1]);
plot([-45 360], [0 0], '-k', 'LineWidth', 0.5);
axis([-45 360 -20 20]); ylabel('Bias (°)'); xlabel('Target Angle (°)');
set(gca, 'xtick', 0:90:360, 'ytick', -20:10:20, 'LineWidth', 1); box on; grid on;
legend('Median \pm SE','Model fit','Location','northwest');
addBIC(BIC_Sober);
%% =========================================================================
%  9. Fit & Plot Sober+TG Model
%  =========================================================================
disp('Fitting Model 7: Sober+TG...');
min_neg_log_likelihood = 1e10;
for i = 1:10
    fun = @(x) motor_biasSober2_VB(data_to_fit, x);
    x0 = [pi/2*rand, pi*2/3*rand, 0.1*rand, 0.1*rand, 8.1*rand+4, 0.01*rand];
    x1 = [0, pi/2, -1, -1, 4, 0]; x2 = [pi/2, pi, 1, 1, inf, inf];
    [x, fval] = fminsearchbnd(fun, x0, x1, x2, options);
    if fval < min_neg_log_likelihood
        min_neg_log_likelihood = fval; best_fit_params_sober_tg = x;
    end
end
BIC_Sober_TG = 2 * min_neg_log_likelihood * 10 + 6 * log(n_data_points * 10);
x = best_fit_params_sober_tg;
l1=3; l2=3; a = x(1); b = x(2); biasa = x(3); biasb = x(4);
x_start_unbiased = cos(a)*l1 + cos(b)*l2; y_start_unbiased = sin(a)*l1 + sin(b)*l2;
x_start_biased = cos(a+biasa)*l1 + cos(b+biasb)*l2; y_start_biased = sin(a+biasa)*l1 + sin(b+biasb)*l2;
[tx, ty] = pol2cart(Ttheta24, 1);
tx = tx + x_start_biased; ty = ty + y_start_biased;
motor_error = zeros(1, 24);
for ta = 1:24
    dist = sqrt(sum([tx(ta)-x_start_unbiased, ty(ta)-y_start_unbiased].^2));
    x_target = x_start_unbiased + (tx(ta)-x_start_unbiased)/dist;
    y_target = y_start_unbiased + (ty(ta)-y_start_unbiased)/dist;
    d = sqrt(x_target^2 + y_target^2);
    theta1 = acos((l1^2 - l2^2 + d^2) / (2*l1*d));
    theta2 = atan2(y_target, x_target);
    a2 = theta2 - theta1; da = a2 - a;
    theta3 = acos((l1^2 + l2^2 - d^2) / (2*l1*l2));
    b2 = pi - theta3 + a2; db = b2 - b;
    rx = cos(a+da+biasa)*l1 + cos(b+db+biasb)*l2;
    ry = sin(a+da+biasa)*l1 + sin(b+db+biasb)*l2;
    [etheta, ~] = cart2pol(rx - x_start_biased, ry - y_start_biased);
    motor_error(ta) = etheta - Ttheta24(ta);
end
p1 = 0:x(5); p2 = (x(5)-0.01):(-2*x(5)/(91-2*x(5)-1)):(-x(5)); p3 = (-x(5)):(-0.1);
target_gen_bias = [p1, p2, p3] * x(6); target_gen_bias = repmat(target_gen_bias, 1, 4);
tg_24_targets = target_gen_bias(1:15:360);
error_theta = motor_error + tg_24_targets;
error_theta(error_theta < -pi) = error_theta(error_theta < -pi) + 2*pi;
figure(107); clf; hold on; title('Model 7: Sober+TG');
errorbar(x24_deg, mdata_deg, se_fit_deg, 'o', ...
    'MarkerFaceColor',[1,0.5,1], 'Color',[1,0.5,1], 'CapSize',0, 'LineStyle','none', 'LineWidth', 1);
plot(x24_deg, error_theta/pi*180, '-', 'linewidth', 2, 'color', [0.5,0.5,1]);
plot([-45 360], [0 0], '-k', 'LineWidth', 0.5);
axis([-45 360 -20 20]); ylabel('Bias (°)'); xlabel('Target Angle (°)');
set(gca, 'xtick', 0:90:360, 'ytick', -20:10:20, 'LineWidth', 1); box on; grid on;
legend('Median \pm SE','Model fit','Location','northwest');
addBIC(BIC_Sober_TG);
