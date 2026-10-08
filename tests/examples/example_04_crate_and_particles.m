%% Example 4 -- C-rate and number of particle-size groups
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

base = default_parameters();
base = set_parameters(base, 'gamma', 1, 'verbose', false);

% C-rate sweep
[sols_c, labels_c] = run_parameter_sweep('c_rate', [0.002, 0.01, 0.05], base, 'n_groups', 10);
plot_voltage_soc(sols_c, labels_c);

% Number of size groups (resolution of the size distribution)
[sols_n, labels_n] = run_parameter_sweep('n_groups', [1, 10, 50], base, 'c_rate', 0.01);
plot_voltage_soc(sols_n, labels_n);
plot_phase_fraction(sols_n, labels_n);

% Different C-rates per segment: slow lithiation, fast delithiation
sol = run_simulation(base, 'n_groups', 10, 'c_rate', [0.005, 0.05]);
plot_voltage_soc(sol);
