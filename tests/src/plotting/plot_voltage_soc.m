function ax = plot_voltage_soc(sols, labels, ax)
%PLOT_VOLTAGE_SOC Electrode potential against state of charge.
%   PLOT_VOLTAGE_SOC(sol) or PLOT_VOLTAGE_SOC({sol1, sol2}, {'a', 'b'})
if isstruct(sols), sols = {sols}; end
if nargin < 2, labels = {}; end
if nargin < 3 || isempty(ax), ax = new_figure(); end
for i = 1:numel(sols)
    plot(ax, sols{i}.SOC, sols{i}.Phi, 'LineWidth', 1.5);
end
xlabel(ax, 'SOC');
ylabel(ax, '$\Phi$ [V]', 'Interpreter', 'latex');
if ~isempty(labels)
    legend(ax, labels, 'Location', 'best', 'Interpreter', 'none');
end
end
