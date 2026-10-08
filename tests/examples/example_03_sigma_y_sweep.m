%% Example 3 -- effect of the yield stress sigma_y of the external shell
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

sigma_ys = [0, 0.05, 0.1, 0.2];
[sols, labels] = run_parameter_sweep('sigma_y', sigma_ys, ...
    'gamma', 1, 'n_groups', 10, 'c_rate', 0.01, 'verbose', false);

plot_voltage_soc(sols, labels);
plot_external_stress(sols{end});