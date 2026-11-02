%% Example 7 -- ensemble with and without external stress (paper Figs. 4-7)
% Ensemble of n = 50 groups at C/500 with mechanics (Lambda = 1.32), without
% external stress and with a rapidly yielding plastic shell, sigma_y = 0.8.
% The shell opposes each particle's surface motion: it tends to -sigma_y
% while the particle expands and +sigma_y while it contracts. After the
% reversal the stress unloads elastically over a change in filling
% fraction Delta c = 2*sigma_y/(Lambda*(1+3*nu_bar)), about 0.48 here, so
% the full gap 2*sigma_y*R_g*T/F is only reached on the second half of
% delithiation.
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
p = sol_with_ext.params;
fprintf(['\nSwitching window Delta c = %.3f. Gap increase at SoC = 0.25, 0.5, 0.75: %s mV ' ...
         '(full value 2*sigma_y*R_g*T/F = %.2f mV, reached outside the window)\n'], ...
    2*p.sigma_y/(p.Lambda*(1 + 3*p.nu_bar)), mat2str(round(s(2).gap - s(1).gap, 2)), ...
    2*p.sigma_y*1000*p.R_gT/p.F);
plot_external_stress(sol_with_ext);