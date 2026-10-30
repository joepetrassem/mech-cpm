function ax = plot_voltage_time(sols, labels, ax)
%PLOT_VOLTAGE_TIME Electrode potential Phi [mV] against time [h].
%   PLOT_VOLTAGE_TIME(sol) or PLOT_VOLTAGE_TIME({sol1, sol2}, {'a', 'b'})
%   Dotted vertical lines mark the protocol blocks of the first run.
%
%   See also PLOT_VOLTAGE_SOC, VOLTAGE_SUMMARY.
if isstruct(sols), sols = {sols}; end
if nargin < 2, labels = {}; end
if nargin < 3 || isempty(ax), ax = new_figure(); end
for i = 1:numel(sols)
    plot(ax, sols{i}.t/3600, 1000*sols{i}.Phi, 'LineWidth', 1.5);
end
t_split = sols{1}.protocol.t_split(2:end-1)/3600;
for k = 1:numel(t_split)
    xline(ax, t_split(k), ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
end
xlabel(ax, 'Time [h]');
ylabel(ax, '$\Phi$ [mV]', 'Interpreter', 'latex');
title(ax, run_title(sols{1}.params, numel(sols)));
if ~isempty(labels)
    legend(ax, labels, 'Location', 'best', 'Interpreter', 'none');
end
end

function s = run_title(p, n_runs)
if n_runs == 1
    s = sprintf('\\Lambda = %g, \\sigma_y = %g, n = %d, C/%g', ...
        p.Lambda, p.sigma_y, p.n_groups, 1/p.c_rate(1));
else
    s = 'Electrode potential';
end
end
