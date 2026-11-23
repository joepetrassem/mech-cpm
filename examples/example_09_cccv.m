%% Example 9 -- constant-current constant-voltage (CCCV) cycling
% A CCCV protocol runs each charge or discharge at constant current until the
% electrode potential reaches a cut-off, then holds the potential there while
% the current decays to I_cutoff. During the hold every particle sees the same
% potential, the reaction currents follow directly from Butler-Volmer, and the
% total current is an output.
%
% The defaults mimic a typical silicon half-cell: lithiation at constant
% current to a cut-off just past the plateau (V_cutoff(1) = -50 mV), then a
% CV hold until |I| < C/40; delithiation at constant current to a cut-off well
% past the plateau (V_cutoff(2) = +70 mV) with no hold. Potentials are
% relative to the reference of mu_chem (0 mV is the plateau of the paper's
% mu_chem), not cell voltages.
%
% With CCCV the SoC waypoints only act as backstops: set them close to 0 and 1
% so that the cut-offs end each step.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

common = {'Lambda', 1, 'sigma_y', 0.6, 'n_groups', 20, 'c_rate', 0.2, 'verbose', false};

% Constant current only (between SoC 0.01 and 0.99), for comparison
sols{1} = run_simulation(common{:});
% CCCV with the defaults: CV hold after lithiation only
sols{2} = run_simulation(common{:}, 'protocol_type', 'CCCV', ...
    'soc_waypoints', [0.01, 0.999, 0.001]);
% CCCV with CV holds after both lithiation and delithiation
sols{3} = run_simulation(common{:}, 'protocol_type', 'CCCV', ...
    'soc_waypoints', [0.01, 0.999, 0.001], 'cv_hold', 'both');
labels = {'CC, C/5', 'CCCV, hold after lithiation', 'CCCV, hold after both'};

plot_voltage_time(sols, labels);
plot_current_time(sols, labels);
plot_voltage_soc(sols, labels);
plot_external_stress(sols{3});
voltage_summary(sols, labels);

% Each constant-current (CC) and constant-voltage (CV) step
for i = 2:3
    sol = sols{i};
    fprintf('\n%s (cut-offs %+g / %+g mV, I_cutoff = C/%g)\n', labels{i}, ...
        sol.params.V_cutoff, 1/sol.params.I_cutoff);
    fprintf('  %-4s %8s %10s %10s %11s %11s %10s\n', 'step', 'segment', 't0 [h]', ...
        'dt [h]', 'SoC start', 'SoC end', 'I end [C]');
    for st = sol.protocol.steps
        in = sol.t >= st.t_start & sol.t <= st.t_end;
        k  = find(in);
        fprintf('  %-4s %8d %10.3f %10.3f %11.4f %11.4f %10.4f\n', st.mode, st.segment, ...
            st.t_start/3600, (st.t_end - st.t_start)/3600, sol.SOC(k(1)), sol.SOC(k(end)), ...
            sol.I(k(end))*3600/sol.params.Q);
    end
end
