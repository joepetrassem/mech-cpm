function ax = plot_ocv_voltage(ocv, ax)
%PLOT_OCV_VOLTAGE Quasi-static electrode potential of a single particle.
%
%   PLOT_OCV_VOLTAGE(ocv) plots the equilibrium electrode potential
%   Phi = -mu*R_g*T/F [mV] (zero overpotential) against c_avg for the
%   lithiation (solid) and delithiation (dashed) branches computed by
%   COMPUTE_OCV_HYSTERESIS, one colour per external stress.
%
%   See also COMPUTE_OCV_HYSTERESIS, PLOT_OCV_HYSTERESIS.
if nargin < 2 || isempty(ax), ax = new_figure(); end
cols = lines(numel(ocv.sig_ext));
for j = 1:numel(ocv.sig_ext)
    lbl = sprintf('\\sigma^{ext} = %g', ocv.sig_ext(j));
    plot(ax, ocv.c_avg, ocv.Phi_lith(j, :),   '-',  'Color', cols(j, :), ...
        'LineWidth', 1.5, 'DisplayName', [lbl ', lithiation']);
    plot(ax, ocv.c_avg, ocv.Phi_delith(j, :), '--', 'Color', cols(j, :), ...
        'LineWidth', 1.5, 'DisplayName', [lbl ', delithiation']);
end
xlim(ax, [0 1]);
xlabel(ax, 'c_{avg}', 'Interpreter', 'tex');
ylabel(ax, '$\Phi_{eq}$ [mV]', 'Interpreter', 'latex');
title(ax, sprintf('Equilibrium potential, \\Lambda = %g', ocv.params.Lambda));
legend(ax, 'Location', 'best');
end
