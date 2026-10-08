function ax = plot_phase_fraction(sols, labels, ax)
%PLOT_PHASE_FRACTION Number fraction of core-shell particles against SOC.
if isstruct(sols), sols = {sols}; end
if nargin < 2, labels = {}; end
if nargin < 3 || isempty(ax), ax = new_figure(); end
for i = 1:numel(sols)
    N = sols{i}.params.N_part;
    frac = double(sols{i}.s ~= 0)*N/sum(N);
    plot(ax, sols{i}.SOC, frac, 'LineWidth', 1.5);
end
xlabel(ax, 'SOC');
ylabel(ax, 'Fraction of phase-separated particles');
ylim(ax, [0, 1.05]);
if ~isempty(labels)
    legend(ax, labels, 'Location', 'best', 'Interpreter', 'none');
end
end
