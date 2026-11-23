function ax = plot_current_time(sols, labels, ax)
%PLOT_CURRENT_TIME Applied current, as a C-rate, against time [h].
%   PLOT_CURRENT_TIME(sol) or PLOT_CURRENT_TIME({sol1, sol2}, {'a', 'b'})
%   plots I*3600/Q (negative while lithiating). Under a CCCV protocol the
%   constant-voltage holds show as the decaying tails of the current; dashed
%   vertical lines mark where the holds of the first run begin.
%
%   See also PLOT_VOLTAGE_TIME, RUN_SIMULATION.
if isstruct(sols), sols = {sols}; end
if nargin < 2, labels = {}; end
if nargin < 3 || isempty(ax), ax = new_figure(); end
for i = 1:numel(sols)
    plot(ax, sols{i}.t/3600, sols{i}.I*3600/sols{i}.params.Q, 'LineWidth', 1.5);
end
mark_cv_starts(ax, sols{1});
xlabel(ax, 'Time [h]');
ylabel(ax, 'Current [C-rate]');
if ~isempty(labels)
    legend(ax, labels, 'Location', 'best', 'Interpreter', 'none');
end
end

function mark_cv_starts(ax, sol)
if ~isfield(sol.protocol, 'steps'), return; end
for st = sol.protocol.steps
    if strcmp(st.mode, 'CV')
        xline(ax, st.t_start/3600, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    end
end
end
