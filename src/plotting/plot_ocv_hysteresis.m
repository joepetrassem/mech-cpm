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
%     - arrowheads tangent to the path, pointing in the direction of charge.
%   Colours: red for compressive (sig_ext < 0), black for sig_ext = 0, blue
%   for tensile (sig_ext > 0).
%
%   Without phase separation (Lambda >= Lambda_c) the two paths coincide on
%   mu_chem(c_avg) - sig_ext; one arrowhead per path is drawn, on the lower
%   half of c_avg for lithiation and the upper half for delithiation.
%
%   PLOT_OCV_HYSTERESIS(ocv, ax, 'Phi') plots the equilibrium electrode
%   potential Phi = -mu*R_g*T/F [mV] instead.
%
%   Arrowheads are drawn last, at the final axis limits and size, so that
%   they look tangent on screen. The axis limits are then fixed. If the
%   figure is resized or the limits are changed afterwards, call
%   PLOT_OCV_HYSTERESIS again to redraw them.
%
%   See also COMPUTE_OCV_HYSTERESIS, PLOT_GIBBS_MINIMISERS.
if nargin < 2 || isempty(ax), [~, ax] = make_fixed_figure(); end
if nargin < 3 || isempty(quantity), quantity = 'mu'; end
quantity = validatestring(quantity, {'mu', 'Phi'}, mfilename, 'quantity');
hold(ax, 'on'); grid(ax, 'on');

VT = 1000*ocv.params.R_gT/ocv.params.F;
if strcmp(quantity, 'Phi'), scale = -VT; else, scale = 1; end

%% Curves
unstable = ocv.single_state == 2;
n_sig = numel(ocv.sig_ext);
[y_l, y_d] = deal(cell(1, n_sig));
cols = zeros(n_sig, 3);
for j = 1:n_sig
    cols(j, :) = stress_colour(ocv.sig_ext(j));
    y_l{j} = scale*ocv.lith.mu(j, :);
    y_d{j} = scale*ocv.delith.mu(j, :);

    plot(ax, ocv.lith.c,   y_l{j}, '-', 'Color', cols(j, :), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\sigma^{ext} = %g', ocv.sig_ext(j)));
    plot(ax, ocv.delith.c, y_d{j}, '-', 'Color', cols(j, :), 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    plot(ax, ocv.c_avg(unstable), scale*ocv.mu_single(j, unstable), ':', ...
        'Color', cols(j, :), 'LineWidth', 1.2, 'HandleVisibility', 'off');
end

%% Axes (before the arrows, which depend on the final limits and size)
xlim(ax, [0 1]);
xticks(ax, 0:0.2:1);
xlabel(ax, 'c_{avg}', 'Interpreter', 'tex');
if strcmp(quantity, 'mu')
    ylim(ax, [-1.5, 1.5]);
    yticks(ax, -1.5:0.5:1.5);
    ylabel(ax, '\mu', 'Interpreter', 'tex');
else
    ylim(ax, ylim(ax));                      % freeze the automatic limits
    ylabel(ax, '$\Phi_{eq}$ [mV]', 'Interpreter', 'latex');
end
if ocv.phase_separation
    title(ax, sprintf('\\Lambda = %g', ocv.params.Lambda));
else
    title(ax, sprintf('\\Lambda = %g (no phase separation)', ocv.params.Lambda));
end
if n_sig > 1, legend(ax, 'Location', 'best'); end

%% Arrowheads: on the single-phase section before the jump, and on the
%  phase-separated section travelled by that branch only
delete(findall(ax, 'Tag', 'lith_arrow', '-or', 'Tag', 'delith_arrow'));   % findall: arrows are hidden from the legend
drawnow;
L = ocv.lith;  D = ocv.delith;  c_sp = ocv.c_sp;
if ocv.phase_separation
    lith_sections   = {1:find(L.s ~= 0, 1) - 1, find(L.s ~= 0 & L.c >= c_sp(2))};
    delith_sections = {1:find(D.s ~= 0, 1) - 1, find(D.s ~= 0 & D.c <= c_sp(1))};
else                                         % the paths coincide: split them
    lith_sections   = {find(L.c <= 0.5)};
    delith_sections = {find(D.c >= 0.5)};
end
for j = 1:n_sig
    for k = 1:numel(lith_sections)
        draw_tangent_arrow(ax, L.c, y_l{j}, lith_sections{k},   cols(j, :), 'lith_arrow');
        draw_tangent_arrow(ax, D.c, y_d{j}, delith_sections{k}, cols(j, :), 'delith_arrow');
    end
end
end

% -------------------------------------------------------------------------
function draw_tangent_arrow(ax, c, y, idx, col, tag)
% Filled arrowhead at the on-screen midpoint of the path section idx,
% aligned with the path tangent in screen (pixel) coordinates and pointing
% in the direction of increasing index, i.e. the direction of travel.
L_px = 9;  W_px = 4.5;                       % arrowhead length and half-width [px]

xl = xlim(ax);  yl = ylim(ax);
pos = getpixelposition(ax);
sx = pos(3)/diff(xl);                         % pixels per data unit
sy = pos(4)/diff(yl);

idx = idx(y(idx) >= yl(1) & y(idx) <= yl(2));  % visible part only
if numel(idx) < 3, return; end
X = c(idx)*sx;  Y = y(idx)*sy;
arc = [0, cumsum(hypot(diff(X), diff(Y)))];
[~, k] = min(abs(arc - arc(end)/2));
k = min(max(k, 2), numel(idx) - 1);

u = [X(k+1) - X(k-1), Y(k+1) - Y(k-1)];
if norm(u) == 0, return; end
u = u/norm(u);                               % unit tangent on screen
n = [-u(2), u(1)];                           % unit normal on screen
P = [L_px*u; -0.5*L_px*u + W_px*n; -0.5*L_px*u - W_px*n];   % tip first

patch(ax, c(idx(k)) + P(:, 1)/sx, y(idx(k)) + P(:, 2)/sy, col, ...
    'EdgeColor', col, 'HandleVisibility', 'off', 'Tag', tag);
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
