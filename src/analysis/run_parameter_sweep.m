function [sols, labels] = run_parameter_sweep(name, values, varargin)
%RUN_PARAMETER_SWEEP Run RUN_SIMULATION for several values of one parameter.
%
%   [sols, labels] = RUN_PARAMETER_SWEEP('gamma', [0 1 2], 'n_groups', 20)
%   [sols, labels] = RUN_PARAMETER_SWEEP('c_rate', {0.01, [0.01 0.05]}, params)
%
%   values may be a numeric vector (one run per element) or a cell array
%   (one run per cell, e.g. for vector-valued parameters). Remaining
%   arguments are passed to RUN_SIMULATION. labels are suitable for legends.

if isnumeric(values)
    values = num2cell(values);
end
sols   = cell(size(values));
labels = cell(size(values));

for i = 1:numel(values)
    label = sprintf('%s = %s', name, mat2str(values{i}, 4));
    fprintf('=== %s  (%d/%d) ===\n', label, i, numel(values));
    sols{i}   = run_simulation(varargin{:}, name, values{i});
    labels{i} = label;
end
end
