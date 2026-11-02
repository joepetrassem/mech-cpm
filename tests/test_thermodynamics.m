function tests = test_thermodynamics
%TEST_THERMODYNAMICS Consistency checks for the free energy, critical points
%   and the surface chemo-mechanical potential in the paper's scaling
%   (Lambda, stresses in units of R_g*T/V_c).
%   Run with:  results = runtests('tests');
tests = functiontests(localfunctions);
end

function setupOnce(~)
addpath(fileparts(fileparts(mfilename('fullpath'))));
setup_paths();
end

function test_mu_is_derivative_of_f(testCase)
c = linspace(0.02, 0.98, 50);
h = 1e-6;
fd = (f_chem(c + h) - f_chem(c - h))/(2*h);
verifyEqual(testCase, mu_chem(c), fd, 'AbsTol', 1e-5);
end

function test_mu_prime_is_derivative_of_mu(testCase)
c = linspace(0.02, 0.98, 50);
h = 1e-6;
fd = (mu_chem(c + h) - mu_chem(c - h))/(2*h);
verifyEqual(testCase, mu_chem_prime(c), fd, 'RelTol', 1e-5);
end

function test_free_energy_matches_paper(testCase)
% mu_chem(c) = log(c/(1-c)) + 3(1-2c)
c = linspace(0.05, 0.95, 20);
verifyEqual(testCase, mu_chem(c), log(c./(1 - c)) + 3*(1 - 2*c), 'AbsTol', 1e-12);
end

function test_spinodal_points_without_mechanics(testCase)
params = set_parameters(default_parameters(), 'Lambda', 0);
c_sp = find_spinodal_points(params);
verifyEqual(testCase, mu_chem_prime(c_sp), [0 0], 'AbsTol', 1e-8);
end

function test_spinodal_points_with_mechanics(testCase)
% Coherent spinodal: mu_chem'(c) + 2*A_int = 0
params = set_parameters(default_parameters(), 'Lambda', 0.99);
c_sp = find_spinodal_points(params);
verifyEqual(testCase, mu_chem_prime(c_sp) + 2*params.A_int, [0 0], 'AbsTol', 1e-8);
end

function test_critical_coupling_Lambda_c(testCase)
% For mu_chem of the paper, spinodal points exist only for Lambda < Lambda_c = 1.8.
p = default_parameters();
below = find_spinodal_points(set_parameters(p, 'Lambda', 0.99*p.Lambda_c));
verifyTrue(testCase, all(isfinite(below)));
verifyEqual(testCase, mean(below), 0.5, 'AbsTol', 1e-8);    % symmetric mu_chem
above = set_parameters(p, 'Lambda', 1.01*p.Lambda_c);
c_sp = verifyWarning(testCase, @() find_spinodal_points(above), ...
    'find_spinodal_points:stable');
verifyTrue(testCase, all(isnan(c_sp)));
% Above Lambda_c the particle is a solid solution: no error, no warning,
% critical points NaN and phase_separation false.
pc = verifyWarningFree(testCase, @() compute_critical_points(set_parameters(above, 'verbose', false)));
verifyTrue(testCase, all(isnan([pc.c_b, pc.c_sp])));
verifyFalse(testCase, pc.phase_separation);
verifyTrue(testCase, default_parameters().phase_separation);
end

function test_binodal_common_tangent_without_mechanics(testCase)
% With Lambda = 0 the binodal points must satisfy the common-tangent construction.
params = set_parameters(default_parameters(), 'Lambda', 0);
c_b = find_binodal_points(params);
mu = mu_chem(c_b);
verifyEqual(testCase, mu(1), mu(2), 'AbsTol', 1e-4);
slope = (f_chem(c_b(2)) - f_chem(c_b(1)))/diff(c_b);
verifyEqual(testCase, slope, mu(1), 'AbsTol', 1e-4);
end

function test_binodal_methods_agree(testCase)
for Lambda = [0, 0.66, 1.32]
    params = set_parameters(default_parameters(), 'Lambda', Lambda);
    [~, info] = find_binodal_points(params, 'Method', 'both', 'GridPoints', 2000);
    verifyLessThan(testCase, max(abs(info.difference)), 1e-3);
end
end

function test_binodal_independent_of_c_avg(testCase)
% The binodal points do not depend on the average concentration used in
% the minimisation (paper, "Minimisation of the Gibbs free energy").
params = set_parameters(default_parameters(), 'Lambda', 0.99);
a = find_binodal_points(params, 'c_avg', 0.4);
b = find_binodal_points(params, 'c_avg', 0.6);
verifyEqual(testCase, a, b, 'AbsTol', 1e-6);
end

function test_critical_point_ordering(testCase)
params = compute_critical_points(set_parameters(default_parameters(), 'Lambda', 0.33));
verifyLessThan(testCase, params.c_b(1), params.c_sp(1));
verifyLessThan(testCase, params.c_sp(2), params.c_b(2));
end

function test_binodal_method_option_is_used(testCase)
a = compute_critical_points(set_parameters(default_parameters(), 'Lambda', 0.66, ...
    'binodal_method', 'analytical'));
g = compute_critical_points(set_parameters(default_parameters(), 'Lambda', 0.66, ...
    'binodal_method', 'gridsearch'));
verifyEqual(testCase, a.c_b, g.c_b, 'AbsTol', 2e-3);
end

function test_mu_hys_matches_paper_equations(testCase)
% mu_single = mu_chem(c_avg) - sig_ext
% mu_sep    = mu_chem(c_s) - 2*A_int*(c_avg - c_s) - sig_ext
p = compute_critical_points(set_parameters(default_parameters(), 'Lambda', 0.99));
sig = 0.6;  c_avg = 0.4;

mu0 = mu_hys(c_avg, 0, sig, p)/p.R_gT;
verifyEqual(testCase, mu0, mu_chem(c_avg) - sig, 'AbsTol', 1e-12);

c_s = p.c_b(2);                                        % lithium-rich shell, s = +1
mu1 = mu_hys(c_avg, 1, sig, p)/p.R_gT;
verifyEqual(testCase, mu1, mu_chem(c_s) - 2*p.A_int*(c_avg - c_s) - sig, 'AbsTol', 1e-12);

c_s = p.c_b(1);                                        % lithium-poor shell, s = -1
mum = mu_hys(c_avg, -1, sig, p)/p.R_gT;
verifyEqual(testCase, mum, mu_chem(c_s) - 2*p.A_int*(c_avg - c_s) - sig, 'AbsTol', 1e-12);
end

function test_mu_sep_equals_core_and_shell_potentials(testCase)
% At equilibrium mu_core = mu_shell = lambda - sig_ext (paper, Appendix E).
p = compute_critical_points(set_parameters(default_parameters(), 'Lambda', 0.99));
sig = -0.3;  c_c = p.c_b(1);  c_s = p.c_b(2);
for c_avg = [0.3, 0.5, 0.7]
    psi = (c_avg - c_s)/(c_c - c_s);
    mu_core  = mu_chem(c_c) - 2*p.A_int*(c_s - c_c)*(1 - psi) - sig;
    mu_shell = mu_chem(c_s) + 2*p.A_int*(c_s - c_c)*psi - sig;
    lambda   = psi*mu_chem(c_c) + (1 - psi)*mu_chem(c_s);
    verifyEqual(testCase, mu_core, mu_shell, 'AbsTol', 1e-6);
    verifyEqual(testCase, mu_shell, lambda - sig, 'AbsTol', 1e-6);
    verifyEqual(testCase, mu_hys(c_avg, 1, sig, p)/p.R_gT, mu_shell, 'AbsTol', 1e-12);
end
end

function test_external_stress_is_rigid_unit_shift(testCase)
% The external stress translates the potential by -sig_ext, independently
% of Lambda and of the phase state (unit coefficient in this scaling).
c_avg = linspace(0.1, 0.9, 9);
for Lambda = [0, 0.99]
    p = compute_critical_points(set_parameters(default_parameters(), 'Lambda', Lambda));
    for s = [0, 1, -1]
        d = (mu_hys(c_avg, s, 0.6, p) - mu_hys(c_avg, s, 0, p))/p.R_gT;
        verifyEqual(testCase, d, -0.6*ones(size(c_avg)), 'AbsTol', 1e-12);
    end
end
end

function test_ocv_equilibrium_potential(testCase)
% Phi_eq = -mu*R_g*T/F on every branch and path, and the external stress
% shifts it by +sig_ext*R_g*T/F.
p = set_parameters(default_parameters(), 'Lambda', 0.99);
ocv = compute_ocv_hysteresis(p, [-0.6, 0, 0.6], 101);
VT = 1000*p.R_gT/p.F;
verifyEqual(testCase, ocv.Phi_lith,      -VT*ocv.mu_lith,      'AbsTol', 1e-12);
verifyEqual(testCase, ocv.Phi_delith,    -VT*ocv.mu_delith,    'AbsTol', 1e-12);
verifyEqual(testCase, ocv.lith.Phi,      -VT*ocv.lith.mu,      'AbsTol', 1e-12);
verifyEqual(testCase, ocv.delith.Phi,    -VT*ocv.delith.mu,    'AbsTol', 1e-12);
verifyEqual(testCase, ocv.Phi_lith(3, :) - ocv.Phi_lith(2, :), 0.6*VT*ones(1, 101), 'AbsTol', 1e-9);
verifyFalse(testCase, isfield(ocv, 'mu_mixed'));
end

function test_ocv_branch_stability(testCase)
% Single-phase states are stable outside [c_b1, c_b2], metastable between a
% binodal and the adjacent spinodal point, and unstable between the spinodals.
ocv = compute_ocv_hysteresis(set_parameters(default_parameters(), 'Lambda', 0.99), 0, 401);
c = ocv.c_avg;  b = ocv.c_b;  sp = ocv.c_sp;
verifyTrue(testCase, all(ocv.single_state(c < b(1) | c > b(2)) == 0));
verifyTrue(testCase, all(ocv.single_state(c > b(1) & c < sp(1)) == 1));
verifyTrue(testCase, all(ocv.single_state(c > sp(2) & c < b(2)) == 1));
verifyTrue(testCase, all(ocv.single_state(c > sp(1) & c < sp(2)) == 2));
% Between the spinodals mu_chem' < -2*A_int <= 0, so the single-phase
% branch decreases there
u = ocv.single_state == 2;
verifyTrue(testCase, all(diff(ocv.mu_single(u)) < 0));
verifyTrue(testCase, all(isnan(ocv.mu_sep(c < b(1) | c > b(2)))));
end

function test_ocv_rich_and_poor_shell_coincide(testCase)
% At equilibrium mu_core = mu_shell, so the rich-shell (s = +1) and
% poor-shell (s = -1) mixed-phase branches are the same line.
for Lambda = [0, 0.99, 1.32]
    p = compute_critical_points(set_parameters(default_parameters(), 'Lambda', Lambda));
    c = linspace(p.c_b(1), p.c_b(2), 50);
    verifyEqual(testCase, mu_hys(c, 1, 0.3, p), mu_hys(c, -1, 0.3, p), 'AbsTol', 1e-6*p.R_gT);
end
end

function test_ocv_cycling_paths(testCase)
% Paper Fig. 3c,d: lithiation jumps down onto the mixed-phase line at c_s1
% and rejoins the single-phase branch continuously at c_b2; delithiation
% jumps up at c_s2 and rejoins continuously at c_b1.
for Lambda = [0, 0.99]
    ocv = compute_ocv_hysteresis(set_parameters(default_parameters(), 'Lambda', Lambda), [-0.6, 0, 0.6]);
    L = ocv.lith;  D = ocv.delith;
    verifyTrue(testCase, all(diff(L.c) >= 0));           % travel order
    verifyTrue(testCase, all(diff(D.c) <= 0));
    verifyEqual(testCase, unique(L.s), [0, 1]);
    verifyEqual(testCase, unique(D.s), [-1, 0]);

    i = find(L.s ~= 0, 1);                               % lithiation jump
    verifyEqual(testCase, L.c(i), ocv.c_sp(1), 'AbsTol', 1e-12);
    verifyEqual(testCase, L.c(i - 1), L.c(i), 'AbsTol', 1e-12);
    verifyLessThan(testCase, L.mu(:, i) - L.mu(:, i - 1), -1e-3);
    e = find(L.s ~= 0, 1, 'last');                       % exit at c_b2
    verifyEqual(testCase, L.mu(:, e), L.mu(:, e + 1), 'AbsTol', 1e-6);

    i = find(D.s ~= 0, 1);                               % delithiation jump
    verifyEqual(testCase, D.c(i), ocv.c_sp(2), 'AbsTol', 1e-12);
    verifyGreaterThan(testCase, D.mu(:, i) - D.mu(:, i - 1), 1e-3);
    e = find(D.s ~= 0, 1, 'last');                       % exit at c_b1
    verifyEqual(testCase, D.mu(:, e), D.mu(:, e + 1), 'AbsTol', 1e-6);

    % Metastable sections lie between the binodal and spinodal points
    verifyTrue(testCase, all(L.c(L.metastable) >= ocv.c_b(1) & L.c(L.metastable) <= ocv.c_sp(1)));
    verifyTrue(testCase, all(D.c(D.metastable) <= ocv.c_b(2) & D.c(D.metastable) >= ocv.c_sp(2)));
end
end

function test_ocv_without_phase_separation(testCase)
% For Lambda >= Lambda_c the particle stays single phase: lithiation and
% delithiation coincide on mu_chem(c_avg) - sig_ext, independently of Lambda,
% every state is stable and there is no mixed-phase line.
ref = [];
for Lambda = [2, 3]          % (exactly at Lambda_c the outcome depends on rounding)
    p = set_parameters(default_parameters(), 'Lambda', Lambda, 'verbose', false);
    ocv = compute_ocv_hysteresis(p, [-0.6, 0, 0.6], 301);
    verifyFalse(testCase, ocv.phase_separation);
    verifyEqual(testCase, ocv.mu_lith, ocv.mu_delith, 'AbsTol', 1e-14);
    verifyEqual(testCase, ocv.mu_lith, mu_chem(ocv.c_avg) - ocv.sig_ext, 'AbsTol', 1e-12);
    verifyTrue(testCase, all(ocv.single_state == 0));
    verifyTrue(testCase, all(isnan(ocv.mu_sep(:))));
    verifyTrue(testCase, all(ocv.lith.s == 0) && all(ocv.delith.s == 0));
    verifyTrue(testCase, all(diff(ocv.lith.c) > 0) && all(diff(ocv.delith.c) < 0));
    verifyEqual(testCase, ocv.lith.mu, fliplr(ocv.delith.mu), 'AbsTol', 1e-12);
    if isempty(ref), ref = ocv.mu_lith; end
    verifyEqual(testCase, ocv.mu_lith, ref, 'AbsTol', 1e-14);     % independent of Lambda
end
end
