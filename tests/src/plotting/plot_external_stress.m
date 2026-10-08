function ax = plot_external_stress(sol, stride, ax)
%PLOT_EXTERNAL_STRESS External radial stress sig_ext of each size group against time.
%   PLOT_EXTERNAL_STRESS(sol, stride) plots every stride-th group (default: ~10 curves).
n = sol.params.n_groups;
if nargin < 2 || isempty(stride), stride = max(1, round(n/10)); end
if nargin < 3 || isempty(ax), ax = new_figure(); end
plot(ax, sol.t/3600, sol.sig_ext(:, 1:stride:end), 'LineWidth', 1);
xlabel(ax, 'Time [h]');
ylabel(ax, '$\sigma^{ext}_{rr}$', 'Interpreter', 'latex');
title(ax, sprintf('\\sigma_y = %g', sol.params.sigma_y));
end
