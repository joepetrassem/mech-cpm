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
verifyError(testCase, @() compute_critical_points(above), 'MechoCPM:noSpinodal');
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
% Phi_eq = -mu*R_g*T/F on every branch, and the external stress shifts it
% by +sig_ext*R_g*T/F.
p = set_parameters(default_parameters(), 'Lambda', 0.99);
ocv = compute_ocv_hysteresis(p, [-0.6, 0, 0.6], 101);
VT = 1000*p.R_gT/p.F;
verifyEqual(testCase, ocv.Phi_lith,   -VT*ocv.mu_lith,   'AbsTol', 1e-12);
verifyEqual(testCase, ocv.Phi_delith, -VT*ocv.mu_delith, 'AbsTol', 1e-12);
verifyEqual(testCase, ocv.Phi_mixed,  -VT*ocv.mu_mixed,  'AbsTol', 1e-12);
verifyEqual(testCase, ocv.Phi_lith(3, :) - ocv.Phi_lith(2, :), 0.6*VT*ones(1, 101), 'AbsTol', 1e-9);
end

