function summary = voltage_summary(sols, labels, soc_points)
%VOLTAGE_SUMMARY Electrode potential on each branch and the hysteresis gap.
%
%   summary = VOLTAGE_SUMMARY(sol) prints, and returns, the electrode
%   potential Phi [mV] at a few states of charge on the lithiation and
%   delithiation branches of a simulation, and the voltage gap between them.
%   VOLTAGE_SUMMARY({sol1, sol2}, {'a', 'b'}) does the same for several runs.
%   VOLTAGE_SUMMARY(..., soc_points) sets the SoC values (default [0.25 0.5 0.75]).
%
%   The first lithiation block (I < 0) and the first delithiation block
%   (I > 0) of the protocol are used. Points outside a branch's SoC range
%   are NaN. Output struct (one element per run):
%     label, soc         run label, SoC values
%     Phi_lith           Phi on the lithiation branch [mV]
%     Phi_delith         Phi on the delithiation branch [mV]
%     gap                Phi_delith - Phi_lith [mV]
%     Phi_range          [min max] of Phi over the whole run [mV]
%
%   With an external stress switching between -sigma_y (lithiation) and
%   +sigma_y (delithiation), the gap grows by 2*sigma_y*R_g*T/F relative to
%   the same run with sigma_y = 0.
%
%   See also PLOT_VOLTAGE_SOC, PLOT_VOLTAGE_TIME.

if isstruct(sols), sols = {sols}; end
if nargin < 2 || isempty(labels)
    labels = arrayfun(@(i) sprintf('run %d', i), 1:numel(sols), 'UniformOutput', false);
end
if nargin < 3 || isempty(soc_points), soc_points = [0.25, 0.5, 0.75]; end
soc_points = soc_points(:)';

summary = struct('label', {}, 'soc', {}, 'Phi_lith', {}, 'Phi_delith', {}, ...
                 'gap', {}, 'Phi_range', {});
for i = 1:numel(sols)
    sol = sols{i};
    r.label      = labels{i};
    r.soc        = soc_points;
    r.Phi_lith   = branch_potential(sol, -1, soc_points);
    r.Phi_delith = branch_potential(sol, +1, soc_points);
    r.gap        = r.Phi_delith - r.Phi_lith;
    r.Phi_range  = 1000*[min(sol.Phi), max(sol.Phi)];
    summary(i)   = r;
end

print_summary(summary);
end

% -------------------------------------------------------------------------
function Phi = branch_potential(sol, direction, soc_points)
% Phi [mV] interpolated at soc_points on the first block with sign(I) = direction.
Phi = NaN(size(soc_points));
I_block = sol.protocol.I;
b = find(sign(I_block) == direction, 1);
if isempty(b), return; end
t_split = sol.protocol.t_split;
in = sol.t >= t_split(b) & sol.t <= t_split(b + 1);
soc = sol.SOC(in);
phi = 1000*sol.Phi(in);
[soc, iu] = unique(soc);
phi = phi(iu);
if numel(soc) < 2, return; end
Phi = interp1(soc, phi, soc_points, 'linear', NaN);
end

function print_summary(summary)
fprintf('\nElectrode potential Phi [mV]\n');
for i = 1:numel(summary)
    r = summary(i);
    fprintf('  %s  (range %.1f to %.1f mV)\n', r.label, r.Phi_range);
    fprintf('    %8s %12s %12s %10s\n', 'SoC', 'lithiation', 'delithiation', 'gap');
    for j = 1:numel(r.soc)
        fprintf('    %8.2f %12.2f %12.2f %10.2f\n', ...
            r.soc(j), r.Phi_lith(j), r.Phi_delith(j), r.gap(j));
    end
end
end
