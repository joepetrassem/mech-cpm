function sol = run_simulation(varargin)
%RUN_SIMULATION Galvanostatic cycling of a phase-changing particle ensemble.
%
%   sol = RUN_SIMULATION()                           default parameters
%   sol = RUN_SIMULATION('Lambda', 1.32, 'sigma_y', 0.8, 'c_rate', 1/500, 'n_groups', 50)
%   sol = RUN_SIMULATION(params)                     full parameter struct
%   sol = RUN_SIMULATION(params, 'Lambda', 0.99)     struct plus overrides
%
%   Any field of DEFAULT_PARAMETERS can be overridden by name (see
%   SET_PARAMETERS for the accepted aliases).
%
%   The model is the ensemble DAE of the paper ("Ensemble of interacting
%   particles"): average
%   concentrations c_avg^(i), reaction currents j_tr^(i), external stresses
%   sig_ext^(i) and the electrode potential Phi, solved with ODE15I. Phase
%   changes are located with event functions (SOLVE_BLOCK).
%
%   The external stress of each particle starts at zero and follows the
%   motion of that particle's surface (SIG_IDEAL, SURFACE_VELOCITY): it is
%   continuous across current blocks, and after a reversal it unloads
%   elastically from one yield value to the other over a change in filling
%   fraction of about 2*sigma_y/(Lambda*(1 + 3*nu_bar)).
%
%   Output fields
%     t        time [s]                             (nt x 1)
%     SOC      state of charge                      (nt x 1)
%     Phi      electrode potential [V]              (nt x 1)
%     I        applied (ramped) current [A]         (nt x 1)
%     c_av     average concentration [mol/m^3]      (nt x n_groups)
%     theta    c_av/c_m                             (nt x n_groups)
%     jtr      reaction current density [A/m^2]     (nt x n_groups)
%     sig_ext  external radial stress, V_c*sigma/(R_g*T) (nt x n_groups)
%     mu       chemo-mechanical potential [J/mol]   (nt x n_groups)
%     s        phase state (0, +1, -1)              (nt x n_groups)
%     params, protocol, n_events, runtime

if nargin >= 1 && isstruct(varargin{1})
    params = varargin{1};
    varargin(1) = [];
else
    params = default_parameters();
end
params = set_parameters(params, varargin{:});   % validates as well

timer    = tic;
params   = build_particle_population(params);
params   = compute_critical_points(params);
protocol = build_current_protocol(params);
n        = params.n_groups;

check_stress_regime(params);

if params.verbose
    fprintf('MechoCPM: %d size groups, Lambda = %g (A_int = %.4g), sigma_y = %g\n', ...
        n, params.Lambda, params.A_int, params.sigma_y);
    if params.phase_separation
        fprintf('  binodal points [%.4f %.4f], spinodal points [%.4f %.4f]\n', ...
            params.c_b, params.c_sp);
    else
        fprintf('  no phase separation (Lambda >= Lambda_c = %.4g): solid solution\n', ...
            params.Lambda_c);
    end
end

% Initial state: uniform concentration, single-phase particles, no external stress
c0      = params.soc_waypoints(1)*params.c_m*ones(n, 1);
s       = zeros(n, 1);
sig_ext = zeros(n, 1);
Phi0    = -mean(mu_hys(c0/params.c_m, s, sig_ext, params))/params.F;
w       = [c0; zeros(n, 1); sig_ext; Phi0];

t_all = []; w_all = []; s_all = []; I_all = [];
n_events = 0;

for b = 1:numel(protocol.I)
    I = protocol.I(b);
    if b == 1, I_prev = 0; else, I_prev = protocol.I(b - 1); end

    if I_prev ~= 0 && sign(I) ~= sign(I_prev)
        params.tau_I = params.tau_I_reversal;
    else
        params.tau_I = params.tau_I_default;
    end

    sig_ext = w(2*n+1:3*n);                  % continuous across blocks

    tspan = protocol.t_split(b:b+1);
    [w0, wdot0] = consistent_initial_conditions(tspan(1), w(1:n), s, ...
        sig_ext, w(end), I, I_prev, tspan(1), params);

    if params.verbose
        res0 = dae_residual(tspan(1), w0, wdot0, s, I, I_prev, tspan(1), params);
        fprintf('  block %d/%d: I = %+.3e A, %.4g h, |res(t0)| = %.1e\n', ...
            b, numel(protocol.I), I, diff(tspan)/3600, norm(res0));
    end

    [t_b, w_b, s_b, ne] = solve_block(tspan, w0, wdot0, s, I, I_prev, params);
    n_events = n_events + ne;

    if params.verbose
        fprintf('           Phi %+.2f -> %+.2f mV (min %+.2f, max %+.2f), %d events\n', ...
            1000*w_b(1, end), 1000*w_b(end, end), 1000*min(w_b(:, end)), ...
            1000*max(w_b(:, end)), ne);
    end

    w = w_b(end, :)';
    s = s_b(end, :)';

    t_all = [t_all; t_b];                                                 %#ok<AGROW>
    w_all = [w_all; w_b];                                                 %#ok<AGROW>
    s_all = [s_all; s_b];                                                 %#ok<AGROW>
    I_all = [I_all; ramp_current(t_b, tspan(1), I_prev, I, params.tau_I)]; %#ok<AGROW>
end

sol = postprocess(t_all, w_all, s_all, I_all, params, protocol);
sol.n_events = n_events;
sol.runtime  = toc(timer);

if params.verbose
    fprintf('  done: %d phase-change events, Phi in [%+.2f, %+.2f] mV, %.1f s\n', ...
        n_events, 1000*min(sol.Phi), 1000*max(sol.Phi), sol.runtime);
end
end

function check_stress_regime(params)
% The plastic-shell model behaves as a rigid-perfectly plastic shell when
% u0 << swelling rate and Lambda*(1+3*nu_bar)*u0 << k*sigma_y. Warn otherwise.
if params.Lambda == 0 || params.sigma_y == 0, return; end
rate = min(params.c_rate)/3600;                    % filling-fraction rate [1/s]
K = params.Lambda*(1 + 3*params.nu_bar);
if params.u0 > 0.1*rate
    warning('MechoCPM:u0Regime', ...
        ['u0 = %g 1/s is not small compared with the swelling rate C-rate/3600 = %g 1/s: ' ...
         'the external stress will not reach +-sigma_y. Reduce u0.'], params.u0, rate);
end
if K*params.u0 > 0.1*params.k*params.sigma_y
    warning('MechoCPM:u0Regime', ...
        ['Lambda*(1+3*nu_bar)*u0 = %g is not small compared with k*sigma_y = %g: ' ...
         'the shell will not behave as rigid-plastic. Reduce u0 or increase k.'], ...
        K*params.u0, params.k*params.sigma_y);
end
end

function sol = postprocess(t, w, s, I, params, protocol)
n       = params.n_groups;
c_av    = w(:, 1:n);
jtr     = w(:, n+1:2*n);
sig_ext = w(:, 2*n+1:3*n);
Phi     = w(:, end);
theta   = c_av/params.c_m;
SOC     = c_av*(params.N_part.*params.F.*params.Vol)/params.Q;

mu = zeros(size(c_av));
for i = 1:numel(t)
    mu(i, :) = mu_hys(theta(i, :), s(i, :), sig_ext(i, :), params);
end

sol = struct('t', t, 'SOC', SOC, 'Phi', Phi, 'I', I, ...
    'c_av', c_av, 'theta', theta, 'jtr', jtr, 'sig_ext', sig_ext, ...
    'mu', mu, 's', s, 'params', params, 'protocol', protocol);
end
