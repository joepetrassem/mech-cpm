function ax = plot_ocv_hysteresis(ocv, ax, show_branches)
%PLOT_OCV_HYSTERESIS Plot output of COMPUTE_OCV_HYSTERESIS (Fig. 3a,b).
%
%   PLOT_OCV_HYSTERESIS(ocv) plots the single-phase curve mu_chem (dark blue)
%   and the mixed-phase equilibrium curve for every external stress
%   (red for sig_ext < 0, blue otherwise).
%   PLOT_OCV_HYSTERESIS(ocv, ax, true) additionally draws the full lithiation
%   (solid) and delithiation (dashed) branches.
if nargin < 2 || isempty(ax), [~, ax] = make_fixed_figure(); hold(ax, 'on'); end
if nargin < 3 || isempty(show_branches), show_branches = false; end
grid(ax, 'on');

plot(ax, ocv.c_avg, ocv.mu_chem, 'Color', [0 0 0.4], 'LineWidth', 1.5);
for j = 1:numel(ocv.sig_ext)
    if ocv.sig_ext(j) < 0, col = 'r'; else, col = 'b'; end
    plot(ax, ocv.c_avg, ocv.mu_mixed(j, :), col, 'LineWidth', 1.5);
    if show_branches
        plot(ax, ocv.c_avg, ocv.mu_lith(j, :),   col, 'LineStyle', '-.',  'LineWidth', 1);
        plot(ax, ocv.c_avg, ocv.mu_delith(j, :), col, 'LineStyle', '--', 'LineWidth', 1);
    end
end
xlim(ax, [0 1]);
ylim(ax, [-1.5, 1.5]);
xticks(ax, 0:0.2:1);
yticks(ax, -1.5:0.5:1.5);
xlabel(ax, 'c_{avg}', 'Interpreter', 'tex');
ylabel(ax, '\mu', 'Interpreter', 'tex');
end
