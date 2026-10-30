function ax = plot_voltage_soc(sols, labels, ax)
%PLOT_VOLTAGE_SOC Electrode potential Phi [mV] against state of charge.
%   PLOT_VOLTAGE_SOC(sol) or PLOT_VOLTAGE_SOC({sol1, sol2}, {'a', 'b'})
%
%   See also PLOT_VOLTAGE_TIME, VOLTAGE_SUMMARY.
if isstruct(sols), sols = {sols}; end
if nargin < 2, labels = {}; end
if nargin < 3 || isempty(ax), ax = new_figure(); end
for i = 1:numel(sols)
    plot(ax, sols{i}.SOC, 1000*sols{i}.Phi, 'LineWidth', 1.5);
end
xlabel(ax, 'SOC');
ylabel(ax, '$\Phi$ [mV]', 'Interpreter', 'latex');
if numel(sols) == 1
    p = sols{1}.params;
    title(ax, sprintf('\\Lambda = %g, \\sigma_y = %g, n = %d', p.Lambda, p.sigma_y, p.n_groups));
end
if ~isempty(labels)
    legend(ax, labels, 'Location', 'best', 'Interpreter', 'none');
end
end
