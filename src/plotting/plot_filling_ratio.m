function [fig, ax] = plot_filling_ratio(sol, highlight, ax)
%PLOT_FILLING_RATIO Filling ratio c_avg^(i)/c_m of every group against time (Figs. 4, 6).
%
%   PLOT_FILLING_RATIO(sol) draws all particle groups in red.
%   PLOT_FILLING_RATIO(sol, [10 40]) additionally highlights groups i = 10
%   (black) and i = 40 (grey) and adds a legend for them; indices beyond
%   n_groups are ignored.
if nargin < 2, highlight = []; end
if nargin < 3 || isempty(ax)
    [fig, ax] = make_fixed_figure();
else
    fig = ancestor(ax, 'figure');
end
hold(ax, 'on'); grid(ax, 'on');

t = sol.t/3600;
plot(ax, t, sol.theta, 'r-', 'LineWidth', 1.0, 'HandleVisibility', 'off');

shades = {[0 0 0], [0.4 0.4 0.4]};
highlight = highlight(highlight >= 1 & highlight <= sol.params.n_groups);
for j = 1:numel(highlight)
    col = shades{min(j, numel(shades))};
    plot(ax, t, sol.theta(:, highlight(j)), 'Color', col, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('i=%d', highlight(j)));
end
if ~isempty(highlight), legend(ax, 'show'); end
xlabel(ax, 'Time (h)');
ylabel(ax, 'c_{avg}^{(i)}', 'Interpreter', 'tex');
end
