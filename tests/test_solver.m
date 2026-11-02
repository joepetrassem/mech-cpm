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
p = set_parameters(default_parameters(), 'alph', 0.4, 'lambda', 0.99, 'sig_y', 0.3);
verifyEqual(testCase, [p.alpha, p.Lambda, p.sigma_y], [0.4, 0.99, 0.3]);
end

function test_old_gamma_name_errors(testCase)
% gamma belongs to the previous scaling; it must not be silently accepted.
verifyError(testCase, @() set_parameters(default_parameters(), 'gamma', 4), ...
    'MechoCPM:set_parameters:renamed');
verifyError(testCase, @() run_simulation('gamm', 4, 'verbose', false), ...
    'MechoCPM:set_parameters:renamed');
end

function test_derived_mechanical_parameters(testCase)
p = set_parameters(default_parameters(), 'Lambda', 1.32);
verifyEqual(testCase, p.A_int, 1.32*(1 + 3*p.nu_bar)/(3*(1 + p.nu_bar)), 'AbsTol', 1e-14);
verifyEqual(testCase, p.Lambda_c, 1.8, 'AbsTol', 1e-12);    % nu_bar = 0.5 (paper)
verifyError(testCase, @() set_parameters(p, 'A_int', 1), 'MechoCPM:set_parameters:derived');
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
p = set_parameters(default_parameters(), 'n_groups', 5, 'Lambda', 0.33, 'sigma_y', 0.1);
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
sol = run_simulation('n_groups', 3, 'c_rate', 0.5, 'Lambda', 0.165, ...
    'soc_waypoints', [0.01, 0.99], 'verbose', false);
verifyEqual(testCase, sol.SOC(end), 0.99, 'AbsTol', 0.02);
% SOC change must equal the charge passed (I < 0 is lithiation)
dSOC_current = -trapz(sol.t, sol.I)/sol.params.Q;
verifyEqual(testCase, sol.SOC(end) - sol.SOC(1), dSOC_current, 'AbsTol', 1e-3);
verifyTrue(testCase, all(isfinite(sol.Phi)));
end

function test_external_stress_follows_current_direction(testCase)
% The plastic shell opposes the motion of the surface: by the end of
% lithiation every particle is at -sigma_y, by the end of delithiation at
% +sigma_y. The stress starts at zero and is continuous across blocks.
sol = run_simulation('n_groups', 3, 'c_rate', 0.5, 'Lambda', 0.33, 'sigma_y', 0.2, ...
    'verbose', false);
n = sol.params.n_groups;
t_split    = sol.protocol.t_split;                          % [rest, lithiation, delithiation]
end_lith   = find(sol.t < t_split(3), 1, 'last');           % end of lithiation block
end_delith = numel(sol.t);                                  % end of delithiation block
verifyEqual(testCase, sol.sig_ext(1, :), zeros(1, n));
% (within 5e-3: particles that pause or give back lithium while others
% react unload partly, then reload)
verifyEqual(testCase, sol.sig_ext(end_lith,   :), -0.2*ones(1, n), 'AbsTol', 5e-3);
verifyEqual(testCase, sol.sig_ext(end_delith, :), +0.2*ones(1, n), 'AbsTol', 5e-3);
verifyLessThan(testCase, max(abs(diff(sol.sig_ext(end_lith:end_lith+1, :)))), 1e-3);   % no reset
end

function test_regression_against_original_code(testCase)
% Compares with results saved from the original code (Clean MCPM:
% main_dae / parameters.m) for four parameter sets. The electrode potential
% must agree to solver tolerance.
%
% The reference was produced with the previous scaling (gamma, stresses in
% units of 2*G, alpha = 0.33). Converting to the paper's scaling,
%     Lambda  = gamma*alpha,   sigma_y(new) = gamma*sigma_y(old),
% must reproduce it exactly: the change of scaling is not a change of model.
%
% Only the references without external stress are used: the original code
% drove the external stress by the global current and reset it at every
% block, which has been replaced by the surface-velocity model.
alpha_ref = 0.33;
data = load(fullfile(fileparts(mfilename('fullpath')), 'data', 'original_code_reference.mat'));
data.ref = data.ref([data.ref.sigma_y] == 0);
verifyNotEmpty(testCase, data.ref);
for i = 1:numel(data.ref)
    r = data.ref(i);
    sol = run_simulation('Lambda', r.gamma*alpha_ref, 'sigma_y', r.gamma*r.sigma_y, ...
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

function test_voltage_summary_and_stress_gap(testCase)
% voltage_summary returns Phi on both branches. A hysteretic external stress
% of +-sigma_y translates the two branches rigidly, so the gap grows by
% exactly 2*sigma_y*R_g*T/F (paper, "Half-cell lithiation and delithiation").
% Lambda and sigma_y are chosen so that the switching window after the
% reversal, 2*sigma_y/(Lambda*(1+3*nu_bar)) = 0.06, ends well before SoC 0.7.
% A single particle is used: in an ensemble, particles that pause or give
% back lithium while others react unload partly, so the shift is not exact.
common = {'n_groups', 1, 'sigma_R', 0, 'c_rate', 0.5, 'Lambda', 1.32, 'verbose', false};
sol0 = run_simulation(common{:}, 'sigma_y', 0);
sol1 = run_simulation(common{:}, 'sigma_y', 0.1);
s = voltage_summary({sol0, sol1}, {'0', '0.1'}, [0.3, 0.5, 0.7]);
verifyTrue(testCase, all(isfinite([s.Phi_lith, s.Phi_delith])));
verifyEqual(testCase, s(1).gap, s(1).Phi_delith - s(1).Phi_lith, 'AbsTol', 1e-12);
verifyGreaterThan(testCase, s(1).gap, 0);                     % delithiation above lithiation
VT = 1000*sol0.params.R_gT/sol0.params.F;                     % [mV]
verifyEqual(testCase, s(2).gap - s(1).gap, 2*0.1*VT*ones(1, 3), 'AbsTol', 0.05);
verifyEqual(testCase, s(1).Phi_range, 1000*[min(sol0.Phi), max(sol0.Phi)], 'AbsTol', 1e-12);
end

function test_cycle_without_phase_separation(testCase)
% Lambda > Lambda_c: the ensemble runs as a solid solution, with no
% phase-change events, every particle single phase, and charge conserved.
% For a single particle, a switched external stress widens the gap by
% exactly 2*sigma_y*R_g*T/F.
sol3 = run_simulation('n_groups', 3, 'c_rate', 0.5, 'Lambda', 2, 'sigma_y', 0.4, 'verbose', false);
verifyFalse(testCase, sol3.params.phase_separation);
verifyEqual(testCase, sol3.n_events, 0);
verifyTrue(testCase, all(sol3.s(:) == 0));
verifyTrue(testCase, all(isfinite(sol3.Phi)));
verifyEqual(testCase, sol3.SOC(end) - sol3.SOC(1), -trapz(sol3.t, sol3.I)/sol3.params.Q, 'AbsTol', 1e-3);

common = {'n_groups', 1, 'sigma_R', 0, 'c_rate', 0.5, 'Lambda', 2, 'verbose', false};
sol0 = run_simulation(common{:}, 'sigma_y', 0);
sol1 = run_simulation(common{:}, 'sigma_y', 0.4);
s = voltage_summary({sol0, sol1}, {'0', '0.4'}, [0.3, 0.5, 0.7]);
VT = 1000*sol0.params.R_gT/sol0.params.F;
verifyEqual(testCase, s(2).gap - s(1).gap, 2*0.4*VT*ones(1, 3), 'AbsTol', 0.05);
end
s = voltage_summary({sol0, sol1}, {'0', '0.4'}, [0.3, 0.5, 0.7]);
VT = 1000*sol0.params.R_gT/sol0.params.F;
verifyEqual(testCase, s(2).gap - s(1).gap, 2*0.4*VT*ones(1, 3), 'AbsTol', 0.05);
end

function test_external_stress_switching_window(testCase)
% Paper, "A simple model for the mechanical response of the surrounding
% material": while |sig_ext| < sigma_y the shell holds the surface, so
% dsig/dt = -Lambda*(1+3*nu_bar)*(dc_avg/dt)/c_m, and after the reversal the
% stress unloads from -sigma_y to +sigma_y over
% Delta c = 2*sigma_y/(Lambda*(1+3*nu_bar)).
Lambda = 1.32;  sigma_y = 0.4;
sol = run_simulation('n_groups', 1, 'sigma_R', 0, 'c_rate', 0.5, 'Lambda', Lambda, ...
    'sigma_y', sigma_y, 'verbose', false);
K = Lambda*(1 + 3*sol.params.nu_bar);
delith = sol.t >= sol.protocol.t_split(3);
c   = sol.theta(delith);
sig = sol.sig_ext(delith);
c_rev = c(1);                                           % filling fraction at the reversal
verifyEqual(testCase, sig(1), -sigma_y, 'AbsTol', 1e-3);
% Mid-window: elastic unloading, sig = -sigma_y + K*(c_rev - c)
mid = abs(c - (c_rev - sigma_y/K)) < 0.02;
verifyNotEmpty(testCase, find(mid, 1));
verifyEqual(testCase, sig(mid), -sigma_y + K*(c_rev - c(mid)), 'AbsTol', 0.02);
% Beyond the window: yielded at +sigma_y
after = c < c_rev - 2*sigma_y/K - 0.1;
verifyEqual(testCase, sig(after), sigma_y*ones(nnz(after), 1), 'AbsTol', 2e-3);
end

function test_no_external_stress_without_coupling(testCase)
% Lambda = 0: the particle never loads the shell, so sig_ext = 0 for any sigma_y.
sol = run_simulation('n_groups', 3, 'c_rate', 0.5, 'Lambda', 0, 'sigma_y', 0.8, ...
    'verbose', false);
verifyEqual(testCase, sol.sig_ext, zeros(size(sol.sig_ext)), 'AbsTol', 1e-12);
end

function test_stress_rate_solves_implicit_law(testCase)
% solve_stress_rate returns the unique root of the implicit relaxation law.
p = set_parameters(default_parameters(), 'Lambda', 0.99, 'sigma_y', 0.6);
sig   = [-0.6; -0.3; 0; 0.3; 0.6; 0.1];
c_dot = p.c_m*[1e-4; 1e-4; -1e-4; -1e-4; 0; 1e-12];
sd = solve_stress_rate(sig, c_dot, p);
% The residual is steep near the root, so check that the root is bracketed
% to within a few ulps rather than evaluating the residual itself.
g = @(x) x + p.k*(sig - sig_ideal(surface_velocity(x, c_dot, p), p));
d = 8*eps(max(abs(sd), p.k*p.sigma_y));
verifyTrue(testCase, all(g(sd - d) <= 0 & g(sd + d) >= 0));
verifyTrue(testCase, all(sd >= -p.k*(sig + p.sigma_y) & sd <= -p.k*(sig - p.sigma_y)));
end

function test_jacobian_matches_finite_differences(testCase)
% DAE_JACOBIAN agrees with a finite-difference Jacobian of DAE_RESIDUAL at a
% consistent state with some particles phase-separated and stressed.
p = set_parameters(default_parameters(), 'n_groups', 4, 'Lambda', 0.99, 'sigma_y', 0.3, ...
    'c_rate', 0.5, 'verbose', false);
p = compute_critical_points(build_particle_population(p));
n = p.n_groups;
c = [0.2; 0.4; 0.6; 0.8]*p.c_m;
s = [0; 1; -1; 0];
sig = [-0.3; -0.1; 0.2; 0.05];
I = -1e-4;
[w, wp] = consistent_initial_conditions(50, c, s, sig, 0, I, 0, 0, p);
[J, Jp] = dae_jacobian(50, w, wp, s, I, 0, 0, p);
F = @(ww, pp) dae_residual(50, ww, pp, s, I, 0, 0, p);
scale = [p.c_m*ones(n, 1); 1e-3*ones(n, 1); ones(n, 1); 1];
for j = 1:3*n + 1
    h = 1e-6*scale(j);
    e = zeros(3*n + 1, 1);  e(j) = h;
    fd  = (F(w + e, wp) - F(w - e, wp))/(2*h);
    verifyEqual(testCase, full(J(:, j)), fd, 'RelTol', 1e-4, 'AbsTol', 1e-9*max(1, norm(fd, inf)));
end
% wdot columns: the stress rows switch over a surface velocity of order u0,
% so the steps for dc/dt and dsig/dt are taken well inside that width
K = p.Lambda*(1 + 3*p.nu_bar);
hp = [1e-3*p.u0_used*p.c_m*ones(n, 1); 1e-6*ones(n, 1); 1e-3*p.u0_used*K*ones(n, 1); 1e-6];
for j = 1:3*n + 1
    e = zeros(3*n + 1, 1);  e(j) = hp(j);
    fd  = (F(w, wp + e) - F(w, wp - e))/(2*hp(j));
    verifyEqual(testCase, full(Jp(:, j)), fd, 'RelTol', 1e-4, 'AbsTol', 1e-9*max(1, norm(fd, inf)));
end
end
end

