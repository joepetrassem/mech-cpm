function ax = plot_potential(sol, stride, ax)
%PLOT_POTENTIAL Chemo-mechanical potential of each size group against time.
%   PLOT_POTENTIAL(sol, stride) plots every stride-th group (default: ~10 curves).
n = sol.params.n_groups;
if nargin < 2 || isempty(stride), stride = max(1, round(n/10)); end
if nargin < 3 || isempty(ax), ax = new_figure(); end
plot(ax, sol.t/3600, sol.mu(:, 1:stride:end)/sol.params.R_gT, 'LineWidth', 1);
xlabel(ax, 'Time [h]');
ylabel(ax, '$\mu/R_gT$', 'Interpreter', 'latex');
title(ax, sprintf('\\gamma = %g, \\alpha = %g, \\sigma_y = %g', ...
    sol.params.gamma, sol.params.alpha, sol.params.sigma_y));
end
