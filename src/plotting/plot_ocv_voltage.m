function ax = plot_ocv_voltage(ocv, ax)
%PLOT_OCV_VOLTAGE Quasi-static electrode potential of a single particle.
%
%   PLOT_OCV_VOLTAGE(ocv) plots the equilibrium electrode potential
%   Phi = -mu*R_g*T/F [mV] (zero overpotential) along the lithiation and
%   delithiation paths from COMPUTE_OCV_HYSTERESIS, with the same layout as
%   PLOT_OCV_HYSTERESIS: metastable sections and jumps solid, unstable states
%   dotted, arrows for the direction of charge.
%
%   See also COMPUTE_OCV_HYSTERESIS, PLOT_OCV_HYSTERESIS.
if nargin < 2 || isempty(ax), ax = new_figure(); end
ax = plot_ocv_hysteresis(ocv, ax, 'Phi');
title(ax, sprintf('Equilibrium potential, \\Lambda = %g', ocv.params.Lambda));
end
