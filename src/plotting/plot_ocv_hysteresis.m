function ax = plot_ocv_hysteresis(ocv, ax, quantity)
%PLOT_OCV_HYSTERESIS Single-particle cycling hysteresis (paper Fig. 3c,d).
%
%   PLOT_OCV_HYSTERESIS(ocv) plots, for every external stress in ocv (from
%   COMPUTE_OCV_HYSTERESIS), the chemo-mechanical potential mu along the
%   lithiation and delithiation paths of one particle:
%     - single-phase sections, including the metastable parts between the
%       binodal and spinodal points (solid);
%     - the jump onto the mixed-phase line at the spinodal point and the
%       phase-separated section up to the far binodal point (solid);
%     - the unstable single-phase states between the spinodal points (dotted);
%     - arrows for the direction of charge (> lithiation, < delithiation).
%   Colours: red for compressive (sig_ext < 0), black for sig_ext = 0, blue
%   for tensile (sig_ext > 0).
%
%   PLOT_OCV_HYSTERESIS(ocv, ax, 'Phi') plots the equilibrium electrode
%   potential Phi = -mu*R_g*T/F [mV] instead.
%
%   See also COMPUTE_OCV_HYSTERESIS, PLOT_GIBBS_MINIMISERS, PLOT_OCV_VOLTAGE.
if nargin < 2 || isempty(ax), [~, ax] = make_fixed_figure(); end
if nargin < 3 || isempty(quantity), quantity = 'mu'; end
quantity = validatestring(quantity, {'mu', 'Phi'}, mfilename, 'quantity');
hold(ax, 'on'); grid(ax, 'on');

VT = 1000*ocv.params.R_gT/ocv.params.F;
if strcmp(quantity, 'Phi'), scale = -VT; else, scale = 1; end

unstable = ocv.single_state == 2;
for j = 1:numel(ocv.sig_ext)
    col = stress_colour(ocv.sig_ext(j));
    y_l = scale*ocv.lith.mu(j, :);
    y_d = scale*ocv.delith.mu(j, :);

    plot(ax, ocv.lith.c,   y_l, '-', 'Color', col, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\sigma^{ext} = %g', ocv.sig_ext(j)));
    plot(ax, ocv.delith.c, y_d, '-', 'Color', col, 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    plot(ax, ocv.c_avg(unstable), scale*ocv.mu_single(j, unstable), ':', ...
        'Color', col, 'LineWidth', 1.2, 'HandleVisibility', 'off');

    % Arrows on sections travelled by one branch only
    c_b = ocv.c_b;  c_sp = ocv.c_sp;
    draw_arrows(ax, ocv.lith,   y_l, [mean([ocv.lith.c(1),   c_b(1)]), mean([c_sp(2), c_b(2)])], '>', col);
    draw_arrows(ax, ocv.delith, y_d, [mean([ocv.delith.c(1), c_b(2)]), mean([c_b(1), c_sp(1)])], '<', col);
end

xlim(ax, [0 1]);
xticks(ax, 0:0.2:1);
xlabel(ax, 'c_{avg}', 'Interpreter', 'tex');
if strcmp(quantity, 'mu')
    ylim(ax, [-1.5, 1.5]);
    yticks(ax, -1.5:0.5:1.5);
    ylabel(ax, '\mu', 'Interpreter', 'tex');
else
    ylabel(ax, '$\Phi_{eq}$ [mV]', 'Interpreter', 'latex');
end
title(ax, sprintf('\\Lambda = %g', ocv.params.Lambda));
if numel(ocv.sig_ext) > 1, legend(ax, 'Location', 'best'); end
end

% -------------------------------------------------------------------------
function draw_arrows(ax, path, y, targets, marker, col)
% One arrowhead on the first single-phase section of the path (near
% targets(1)) and one on its phase-separated section (near targets(2)).
first_sep = find(path.s ~= 0, 1);
sections = {1:first_sep - 1, find(path.s ~= 0)};
for k = 1:2
    idx = sections{k};
    [~, i] = min(abs(path.c(idx) - targets(k)));
    plot(ax, path.c(idx(i)), y(idx(i)), marker, 'Color', col, ...
        'MarkerFaceColor', col, 'MarkerSize', 5, 'HandleVisibility', 'off');
end
end

function col = stress_colour(sig)
if sig < 0
    col = [0.85 0 0];       % compressive
elseif sig > 0
    col = [0 0.3 0.9];      % tensile
else
    col = [0 0 0];
end
end
