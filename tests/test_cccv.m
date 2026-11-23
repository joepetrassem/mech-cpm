function tests = test_cccv
%TEST_CCCV Constant-current constant-voltage protocol.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
addpath(fileparts(fileparts(mfilename('fullpath'))));
setup_paths();
args = {'n_groups', 5, 'c_rate', 0.2, 'Lambda', 1, 'sigma_y', 0.6, 'verbose', false, ...
        'protocol_type', 'CCCV', 'soc_waypoints', [0.01, 0.999, 0.001]};
testCase.TestData.both = run_simulation(args{:}, 'cv_hold', 'both');
testCase.TestData.lith = run_simulation(args{:});            % default: hold after lithiation
end

function test_cc_is_default(testCase)
p = default_parameters();
verifyEqual(testCase, p.protocol_type, 'CC');
sol = run_simulation('n_groups', 3, 'c_rate', 0.5, 'verbose', false);
verifyFalse(testCase, any(sol.cv));
verifyEqual(testCase, {sol.protocol.steps.mode}, {'CC', 'CC', 'CC'});   % rest, lith, delith
end

function test_potential_held_at_cutoffs(testCase)
sol = testCase.TestData.both;
V = sol.params.V_cutoff/1000;
verifyEqual(testCase, {sol.protocol.steps.mode}, {'CC', 'CC', 'CV', 'CC', 'CV'});
for st = sol.protocol.steps(strcmp({sol.protocol.steps.mode}, 'CV'))
    in = sol.cv & sol.t >= st.t_start & sol.t <= st.t_end;
    lith = sol.I(find(in, 1)) < 0;
    verifyEqual(testCase, sol.Phi(in), V(2 - lith)*ones(nnz(in), 1), 'AbsTol', 1e-12);
end
% the constant-current steps end at the cut-off and never go past it
verifyGreaterThanOrEqual(testCase, min(sol.Phi), V(1) - 1e-6);
verifyLessThanOrEqual(testCase, max(sol.Phi), V(2) + 1e-6);
end

function test_cv_current_decays_to_cutoff(testCase)
sol = testCase.TestData.both;
I_cut = sol.params.I_cutoff*sol.params.Q/3600;
for st = sol.protocol.steps(strcmp({sol.protocol.steps.mode}, 'CV'))
    in = find(sol.cv & sol.t >= st.t_start & sol.t <= st.t_end);
    I = abs(sol.I(in));
    verifyGreaterThan(testCase, I(1), 2*I_cut);                 % starts near C/5
    verifyEqual(testCase, I(end), I_cut, 'RelTol', 1e-3);         % ends at C/40
    verifyLessThanOrEqual(testCase, max(diff(I)), 1e-3*I_cut);    % decays monotonically
end
end

function test_cccv_conserves_charge(testCase)
for sol = {testCase.TestData.both, testCase.TestData.lith}
    s = sol{1};
    dSOC = -trapz(s.t, s.I)/s.params.Q;
    verifyEqual(testCase, s.SOC(end) - s.SOC(1), dSOC, 'AbsTol', 1e-3);
end
end

function test_cv_hold_adds_capacity(testCase)
% The hold after lithiation fills the electrode beyond the end of the CC step.
sol = testCase.TestData.lith;
st = sol.protocol.steps;
verifyEqual(testCase, {st.mode}, {'CC', 'CC', 'CV', 'CC'});
soc_at = @(t) sol.SOC(find(sol.t <= t, 1, 'last'));
verifyGreaterThan(testCase, soc_at(st(3).t_end) - soc_at(st(3).t_start), 1e-3);
% no hold after delithiation: the run ends at the delithiation cut-off
verifyFalse(testCase, any(sol.cv(sol.t > st(4).t_start)));
verifyEqual(testCase, sol.Phi(end), sol.params.V_cutoff(2)/1000, 'AbsTol', 1e-6);
end

function test_cccv_parameters_validated(testCase)
p = default_parameters();
id = 'MechoCPM:params';
verifyError(testCase, @() set_parameters(p, 'protocol_type', 'CVCC'), id);
verifyError(testCase, @() set_parameters(p, 'cv_hold', 'always'), id);
verifyError(testCase, @() set_parameters(p, 'protocol_type', 'CCCV', 'c_rate', 0.2, ...
    'V_cutoff', [70, -50]), id);
verifyError(testCase, @() set_parameters(p, 'protocol_type', 'CCCV', 'c_rate', 0.01), id); % I_cutoff >= c_rate
verifyError(testCase, @() set_parameters(p, 'protocol_type', 'CCCV', 'c_rate', 0.2, ...
    'solver', 'ode15i'), id);
q = set_parameters(p, 'protocol_type', 'CCCV', 'c_rate', 0.2);
verifyEqual(testCase, q.V_cutoff, [-50, 70]);
verifyEqual(testCase, q.I_cutoff, 0.025);
end
