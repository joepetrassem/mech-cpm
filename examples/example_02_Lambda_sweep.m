%% Example 2 -- effect of the mechanical coupling Lambda
% Lambda reshapes the mixed-phase plateau (through A_int) without widening
% the hysteresis gap. Phase separation is suppressed for Lambda >= 1.8.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

Lambdas = [0, 0.33, 0.66, 1.32];
[sols, labels] = run_parameter_sweep('Lambda', Lambdas, ...
    'n_groups', 10, 'c_rate', 0.01, 'verbose', false);

plot_voltage_soc(sols, labels);
plot_voltage_time(sols, labels);
voltage_summary(sols, labels);
plot_phase_fraction(sols, labels);
