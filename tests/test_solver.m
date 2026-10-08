function tests = test_solver
%TEST_SOLVER Initial-condition consistency, parameter handling, a short run
%   and a regression test against the original (pre-reorganisation) code.
tests = functiontests(localfunctions);
end

function setupOnce(~)
addpath(fileparts(fileparts(mfilename('fullpath'))));
setup_paths();
end

function test_aliases(testCase)
p = set_parameters(default_parameters(), 'alph', 0.4, 'Gamma', 2, 'sig_y', 0.3);
verifyEqual(testCase, [p.alpha, p.gamma, p.sigma_y], [0.4, 2, 0.3]);
end

function test_unknown_parameter_errors(testCase)
verifyError(testCase, @() set_parameters(default_parameters(), 'gama', 1), ...
    'MechoCPM:set_parameters:unknown');
end

function test_derived_parameter_errors(testCase)
verifyError(testCase, @() set_parameters(default_parameters(), 'nu_bar', 1), ...
    'MechoCPM:set_parameters:derived');
end

function test_invalid_value_errors(testCase)
verifyError(testCase, @() set_parameters(default_parameters(), 'c_rate', -1), ...
    ?MException);
end

function test_nu_bar_definition(testCase)
p = default_parameters();
verifyEqual(testCase, p.nu_bar, p.nu/(1 - 2*p.nu), 'AbsTol', 1e-14);
verifyEqual(testCase, p.nu_bar, 0.5, 'AbsTol', 1e-14);      % value used in the paper
end

function test_consistent_initial_conditions(testCase)
p = set_parameters(default_parameters(), 'n_groups', 5, 'gamma', 1, 'sigma_y', 0.1);
p = compute_critical_points(build_particle_population(p));
n = p.n_groups;
c = linspace(0.1, 0.9, n)'*p.c_m;
s = [0; 1; 0; -1; 0];
sig = 0.05*ones(n, 1);
I = -1e-4;
[w0, wdot0] = consistent_initial_conditions(10, c, s, sig, 0, I, 0, 0, p);
res = dae_residual(10, w0, wdot0, s, I, 0, 0, p);
verifyLessThan(testCase, norm(res(n+1:end)), 1e-8);
end

function test_short_cycle_conserves_charge(testCase)
sol = run_simulation('n_groups', 3, 'c_rate', 0.5, 'gamma', 0.5, ...
    'soc_waypoints', [0.01, 0.99], 'verbose', false);
verifyEqual(testCase, sol.SOC(end), 0.99, 'AbsTol', 0.02);
% SOC change must equal the charge passed (I < 0 is lithiation)
dSOC_current = -trapz(sol.t, sol.I)/sol.params.Q;
verifyEqual(testCase, sol.SOC(end) - sol.SOC(1), dSOC_current, 'AbsTol', 1e-3);
verifyTrue(testCase, all(isfinite(sol.Phi)));
end

function test_external_stress_follows_current_direction(testCase)
sol = run_simulation('n_groups', 3, 'c_rate', 0.5, 'gamma', 1, 'sigma_y', 0.2, ...
    'verbose', false);
n = sol.params.n_groups;
t_split    = sol.protocol.t_split;                          % [rest, lithiation, delithiation]
mid_lith   = find(sol.t < t_split(3), 1, 'last');           % end of lithiation block
mid_delith = numel(sol.t);                                  % end of delithiation block
verifyEqual(testCase, sol.sig_ext(mid_lith,   :), -0.2*ones(1, n), 'AbsTol', 1e-3);
verifyEqual(testCase, sol.sig_ext(mid_delith, :), +0.2*ones(1, n), 'AbsTol', 1e-3);
end

function test_regression_against_original_code(testCase)
% Compares with results saved from the original code (Clean MCPM:
% main_dae / parameters.m) for four parameter sets. The electrode potential
% must agree to solver tolerance.
data = load(fullfile(fileparts(mfilename('fullpath')), 'data', 'original_code_reference.mat'));
for i = 1:numel(data.ref)
    r = data.ref(i);
    sol = run_simulation('gamma', r.gamma, 'sigma_y', r.sigma_y, ...
        'n_groups', r.n_groups, 'c_rate', r.c_rate, 'verbose', false);
    verifyEqual(testCase, sol.params.c_b,  r.c_b,  'AbsTol', 1e-8);
    verifyEqual(testCase, sol.params.c_sp, r.c_sp, 'AbsTol', 1e-8);

    [t, iu] = unique(sol.t);
    tg   = linspace(0, min(t(end), r.t(end)), 5000);
    dPhi = interp1(r.t, r.Phi, tg) - interp1(t, sol.Phi(iu), tg);
    dPhi = sort(abs(dPhi(~isnan(dPhi))));
    verifyLessThan(testCase, sqrt(mean(dPhi.^2)), 1e-5, ...
        sprintf('Phi differs from original code (case %d)', i));
    verifyLessThan(testCase, dPhi(ceil(0.99*numel(dPhi))), 1e-4);   % 99th percentile
end
end