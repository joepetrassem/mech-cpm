%% Example 7 -- ensemble with and without external stress (paper Figs. 4-7)
% Ensemble of n = 50 groups at C/500 with mechanics (Lambda = 1.32), without
% external stress and with a plastic external shell that switches between
% -sigma_y (lithiation) and +sigma_y (delithiation), sigma_y = 0.8.
% Takes a few minutes: the cycle lasts about 1000 simulated hours.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

n_groups = 50;
c_rate   = 1/500;

sol_no_ext   = run_simulation('Lambda', 1.32, 'sigma_y', 0,   'n_groups', n_groups, 'c_rate', c_rate);
sol_with_ext = run_simulation('Lambda', 1.32, 'sigma_y', 0.8, 'n_groups', n_groups, 'c_rate', c_rate);

plot_results_comparison(sol_no_ext, sol_with_ext);   % Phi vs SoC, filling ratio, Phi vs time
s = voltage_summary({sol_no_ext, sol_with_ext}, {'sigma_y = 0', 'sigma_y = 0.8'});
fprintf('\nGap increase at SoC = 0.5: %.2f mV (predicted 2*sigma_y*R_g*T/F = %.2f mV)\n', ...
    s(2).gap(2) - s(1).gap(2), 2*0.8*1000*sol_no_ext.params.R_gT/sol_no_ext.params.F);