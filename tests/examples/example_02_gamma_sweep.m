%% Example 2 -- effect of the mechanical coupling gamma
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

gammas = [0, 0.5, 1, 2];
[sols, labels] = run_parameter_sweep('gamma', gammas, ...
    'n_groups', 10, 'c_rate', 0.01, 'verbose', false);

plot_voltage_soc(sols, labels);
plot_phase_fraction(sols, labels);