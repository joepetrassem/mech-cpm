%% Example 3 -- effect of the yield stress sigma_y of the external shell
% sigma_y is scaled by R_g*T/V_c, so a hysteretic external stress of
% +-sigma_y widens the voltage gap by 2*sigma_y*R_g*T/F (about 52 mV per unit).
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

sigma_ys = [0, 0.2, 0.4, 0.8];
[sols, labels] = run_parameter_sweep('sigma_y', sigma_ys, ...
    'Lambda', 0.33, 'n_groups', 10, 'c_rate', 0.01, 'verbose', false);

plot_voltage_soc(sols, labels);
plot_external_stress(sols{end});
