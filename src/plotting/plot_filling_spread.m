function ax = plot_filling_spread(sols, labels, ax)
%PLOT_FILLING_SPREAD Spread of particle filling, std of c_avg across size groups, against SOC.
%   PLOT_FILLING_SPREAD(sol) or PLOT_FILLING_SPREAD({sol1, sol2}, {'a', 'b'})
%   plots std(sol.theta, 0, 2), the sample standard deviation of c_avg over
%   the size groups at each output time (unweighted by particle number).
%   It is zero while all groups fill uniformly and grows once the
%   inter-particle (mosaic) instability sets in.
%
%   See also PLOT_PHASE_FRACTION, PLOT_FILLING_RATIO.
if isstruct(sols), sols = {sols}; end
if nargin < 2, labels = {}; end
if nargin < 3 || isempty(ax), ax = new_figure(); end
for i = 1:numel(sols)
    spread = std(sols{i}.theta, 0, 2);
    plot(ax, sols{i}.SOC, spread, 'LineWidth', 1.5);
end
xlabel(ax, 'SOC');
ylabel(ax, 'Standard deviation of $\bar{c}^{(i)}$', 'Interpreter', 'latex');
ylim(ax, [0, Inf]);
if ~isempty(labels)
    legend(ax, labels, 'Location', 'best', 'Interpreter', 'none');
end
end
