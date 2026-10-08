%% Example 7 -- ensemble with and without external stress (paper Figs. 4-7)
% Ensemble of n = 50 groups at C/500 with mechanics (gamma = 4), without
% external stress and with a plastic external shell that switches between
% -sigma_y (lithiation) and +sigma_y (delithiation), sigma_y = 0.2.
% Takes a few minutes: the cycle lasts about 1000 simulated hours.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

n_groups = 50;
c_rate   = 1/500;

sol_no_ext   = run_simulation('gamma', 4, 'sigma_y', 0,   'n_groups', n_groups, 'c_rate', c_rate);
sol_with_ext = run_simulation('gamma', 4, 'sigma_y', 0.2, 'n_groups', n_groups, 'c_rate', c_rate);

plot_results_comparison(sol_no_ext, sol_with_ext);