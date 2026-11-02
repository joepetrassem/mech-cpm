function ax = plot_gibbs_minimisers(ocv, ax, which_stress)
%PLOT_GIBBS_MINIMISERS Single- and mixed-phase solution branches (paper Fig. 3a,b).
%
%   PLOT_GIBBS_MINIMISERS(ocv) plots the chemo-mechanical potential of the
%   Gibbs free energy minimisers from COMPUTE_OCV_HYSTERESIS: the
%   single-phase branch over the whole range (dark blue) and the mixed-phase
%   branch on [c_b1, c_b2] (blue), for the first external stress in ocv.
%   PLOT_GIBBS_MINIMISERS(ocv, ax, j) uses the j-th external stress.
%
%   Where the mixed-phase solution exists it is the global minimiser; the
%   cycling paths, metastable and unstable states are shown by
%   PLOT_OCV_HYSTERESIS.
if nargin < 2 || isempty(ax), [~, ax] = make_fixed_figure(); end
if nargin < 3 || isempty(which_stress), which_stress = 1; end
hold(ax, 'on'); grid(ax, 'on');
j = which_stress;
plot(ax, ocv.c_avg, ocv.mu_single(j, :), 'Color', [0 0 0.45], 'LineWidth', 1.5, ...
    'DisplayName', 'single phase');
if ocv.phase_separation
    plot(ax, ocv.c_avg, ocv.mu_sep(j, :), 'Color', [0.3 0.55 1], 'LineWidth', 1.5, ...
        'DisplayName', 'mixed phase');
end
xlim(ax, [0 1]); ylim(ax, [-1.5, 1.5]);
xticks(ax, 0:0.2:1); yticks(ax, -1.5:0.5:1.5);
xlabel(ax, 'c_{avg}', 'Interpreter', 'tex');
ylabel(ax, '\mu', 'Interpreter', 'tex');
if ocv.phase_separation
    title(ax, sprintf('\\Lambda = %g', ocv.params.Lambda));
else
    title(ax, sprintf('\\Lambda = %g (no mixed-phase solution)', ocv.params.Lambda));
end
legend(ax, 'Location', 'best');
end
