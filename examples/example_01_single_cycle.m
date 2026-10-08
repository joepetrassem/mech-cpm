%% Example 1 -- a single lithiation/delithiation cycle
% Runs one cycle with mechanics and an external stress switched on, then
% plots the voltage curve, particle potentials, external stress and the
% fraction of phase-separated particles.

clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

%% Parameters
% Every field of default_parameters() can be set here by name.
sol = run_simulation( ...
    'alpha',    0.33, ...   % swelling parameter
    'gamma',    1.0,  ...   % mechanical coupling (0 = no mechanics)
    'sigma_y',  0.1,  ...   % yield stress of the external shell
    'n_groups', 20,   ...   % number of particle-size groups
    'c_rate',   0.01);      % C-rate (scalar or one per segment)

%% Plots
plot_voltage_soc(sol);
plot_potential(sol);
plot_external_stress(sol);
plot_phase_fraction(sol);

%% The same run, built from a parameter struct
% params = default_parameters();
% params = set_parameters(params, 'gamma', 1.0, 'sigma_y', 0.1, ...
%                                 'n_groups', 20, 'c_rate', 0.01);
% sol = run_simulation(params);