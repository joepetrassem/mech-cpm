function ax = plot_phase_fraction(sols, labels, ax)
%PLOT_PHASE_FRACTION Number fraction of core-shell particles against SOC.
%
%   PLOT_PHASE_FRACTION(sol) or PLOT_PHASE_FRACTION({sol1, sol2}, {'a', 'b'})
%   Lithiation is drawn in blue and delithiation in orange
%   (DIRECTION_COLOUR); with several runs, each run has its own shade.
if isstruct(sols), sols = {sols}; end
if nargin < 2, labels = {}; end
if nargin < 3 || isempty(ax), ax = new_figure(); end
n_runs = numel(sols);
for i = 1:n_runs
    N    = sols{i}.params.N_part;
    frac = double(sols{i}.s ~= 0)*N/sum(N);
    seg  = direction_segments(sols{i});
    shown = struct('lithiation', false, 'delithiation', false);
    for k = 1:numel(seg)
        d = seg(k).direction;
        if ~shown.(d)
            if isempty(labels)
                name = d;
            else
                name = sprintf('%s, %s', labels{i}, d);
            end
            vis = 'on';  shown.(d) = true;
        else
            name = '';  vis = 'off';
        end
        plot(ax, sols{i}.SOC(seg(k).idx), frac(seg(k).idx), 'LineWidth', 1.5, ...
            'Color', direction_colour(d, i, n_runs), 'DisplayName', name, ...
            'HandleVisibility', vis, 'Tag', d);
    end
end
xlabel(ax, 'SOC');
ylabel(ax, 'Fraction of phase-separated particles');
ylim(ax, [0, 1.05]);
legend(ax, 'Location', 'best', 'Interpreter', 'none');
end
