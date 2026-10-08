function ax = new_figure(width_cm, height_cm)
%NEW_FIGURE Figure of fixed size with consistent styling; returns its axes.
if nargin < 1, width_cm = 14; end
if nargin < 2, height_cm = 10; end
fig = figure('Color', 'w', 'Units', 'centimeters');
fig.Position(3:4) = [width_cm, height_cm];
ax = axes(fig);
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
set(ax, 'FontSize', 12, 'LineWidth', 0.8);
end
