%% Example 1 -- a single lithiation/delithiation cycle
% Runs one cycle with mechanics and an external stress switched on, then
% plots the electrode potential against SoC and time, the particle
% potentials, the external stress and the fraction of phase-separated
% particles, and prints the electrode potential on each branch.

clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

%% Parameters
% Every field of default_parameters() can be set here by name.
sol = run_simulation( ...
    'Lambda',   0.99, ...   % mechanical coupling (0 = no internal stress; < 1.8)
    'sigma_y',  0.6,  ...   % scaled yield stress of the external shell, V_c*sigma_y/(R_g*T)
    'n_groups', 20,   ...   % number of particle-size groups
    'c_rate',   0.01);      % C-rate (scalar or one per segment)

%% Plots
plot_voltage_soc(sol);
plot_voltage_time(sol);
voltage_summary(sol, {'Lambda = 0.99, sigma_y = 0.6'});
plot_potential(sol);
plot_external_stress(sol);
plot_phase_fraction(sol);

%% The same run, built from a parameter struct
% params = default_parameters();
% params = set_parameters(params, 'Lambda', 0.99, 'sigma_y', 0.6, ...
%                                 'n_groups', 20, 'c_rate', 0.01);
% sol = run_simulation(params);