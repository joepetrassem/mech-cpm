function tests = test_thermodynamics
%TEST_THERMODYNAMICS Consistency checks for the free energy and critical points.
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

function test_free_energy_matches_paper_eq72(testCase)
% mu_chem(c) = log(c/(1-c)) + 3(1-2c)   (Eq. 72)
c = linspace(0.05, 0.95, 20);
verifyEqual(testCase, mu_chem(c), log(c./(1 - c)) + 3*(1 - 2*c), 'AbsTol', 1e-12);
end

function test_spinodal_points_without_mechanics(testCase)
params = set_parameters(default_parameters(), 'gamma', 0);
c_sp = find_spinodal_points(params);
verifyEqual(testCase, mu_chem_prime(c_sp), [0 0], 'AbsTol', 1e-8);
end

function test_spinodal_points_with_mechanics(testCase)
% Eq. (71): mu_chem'(c) + 2*gamma*alpha*(1+3nu_bar)/(3(1+nu_bar)) = 0
params = set_parameters(default_parameters(), 'gamma', 3);
c_sp = find_spinodal_points(params);
K = 2*params.gamma*params.alpha*(1 + 3*params.nu_bar)/(3*(1 + params.nu_bar));
verifyEqual(testCase, mu_chem_prime(c_sp) + K, [0 0], 'AbsTol', 1e-8);
end

function test_binodal_common_tangent_without_mechanics(testCase)
% With gamma = 0 the binodal points must satisfy the common-tangent construction.
params = set_parameters(default_parameters(), 'gamma', 0);
c_b = find_binodal_points(params);
mu = mu_chem(c_b);
verifyEqual(testCase, mu(1), mu(2), 'AbsTol', 1e-4);
slope = (f_chem(c_b(2)) - f_chem(c_b(1)))/diff(c_b);
verifyEqual(testCase, slope, mu(1), 'AbsTol', 1e-4);
end

function test_binodal_methods_agree(testCase)
for gamma = [0, 2, 4]
    params = set_parameters(default_parameters(), 'gamma', gamma);
    [~, info] = find_binodal_points(params, 'Method', 'both', 'GridPoints', 2000);
    verifyLessThan(testCase, max(abs(info.difference)), 1e-3);
end
end

function test_critical_point_ordering(testCase)
params = compute_critical_points(set_parameters(default_parameters(), 'gamma', 1));
verifyLessThan(testCase, params.c_b(1), params.c_sp(1));
verifyLessThan(testCase, params.c_sp(2), params.c_b(2));
end

function test_binodal_method_option_is_used(testCase)
a = compute_critical_points(set_parameters(default_parameters(), 'gamma', 2, ...
    'binodal_method', 'analytical'));
g = compute_critical_points(set_parameters(default_parameters(), 'gamma', 2, ...
    'binodal_method', 'gridsearch'));
verifyEqual(testCase, a.c_b, g.c_b, 'AbsTol', 2e-3);
end

function test_mu_hys_matches_paper_equations(testCase)
% Single phase: mu = mu_chem - gamma*sig_ext          (Eq. 20)
% Core-shell  : mu = mu_chem(c_s) - gamma*(sig_ext - 2(c_s-c_c)*alpha*(1+3nu_bar)/(3(1+nu_bar))*psi)  (Eq. 73)
p = compute_critical_points(set_parameters(default_parameters(), 'gamma', 3));
sig = 0.2;  c_avg = 0.4;
k = (1 + 3*p.nu_bar)/(3*(1 + p.nu_bar));

mu0 = mu_hys(c_avg, 0, sig, p)/p.R_gT;
verifyEqual(testCase, mu0, mu_chem(c_avg) - p.gamma*sig, 'AbsTol', 1e-12);

c_c = p.c_b(1); c_s = p.c_b(2);                       % lithium-rich shell, s = +1
psi = (c_avg - c_s)/(c_c - c_s);
mu1 = mu_hys(c_avg, 1, sig, p)/p.R_gT;
verifyEqual(testCase, mu1, mu_chem(c_s) - p.gamma*(sig - 2*(c_s - c_c)*p.alpha*k*psi), 'AbsTol', 1e-12);

c_c = p.c_b(2); c_s = p.c_b(1);                       % lithium-poor shell, s = -1
psi = (c_avg - c_s)/(c_c - c_s);
mum = mu_hys(c_avg, -1, sig, p)/p.R_gT;
verifyEqual(testCase, mum, mu_chem(c_s) - p.gamma*(sig - 2*(c_s - c_c)*p.alpha*k*psi), 'AbsTol', 1e-12);
end