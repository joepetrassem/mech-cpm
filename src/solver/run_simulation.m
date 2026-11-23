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
%   particles"). It is solved, according to params.solver, either as a
%   reduced ODE in (c_avg, sig_ext) with ODE15S (default; SOLVE_BLOCK_ODE15S)
%   or as the fully implicit DAE with ODE15I (SOLVE_BLOCK); see
%   docs/solver_formulation.md. The variables are the average
%   concentrations c_avg^(i), reaction currents j_tr^(i), external stresses
%   sig_ext^(i) and the electrode potential Phi. Phase changes are located
%   with event functions in either case.
%
%   The external stress of each particle starts at zero and follows the
%   motion of that particle's surface (SIG_IDEAL, SURFACE_VELOCITY): it is
%   continuous across current blocks, and after a reversal it unloads
%   elastically from one yield value to the other over a change in filling
%   fraction of about 2*sigma_y/(Lambda*(1 + 3*nu_bar)).
%
%   Protocol. With protocol_type = 'CC' (default) the current follows
%   soc_waypoints and c_rate: one constant-current segment between each pair
%   of waypoints. With protocol_type = 'CCCV' each segment runs at constant
%   current until Phi reaches the cut-off V_cutoff (lithiation: V_cutoff(1),
%   delithiation: V_cutoff(2), in mV) or the SoC waypoint, whichever comes
%   first; if it stopped at the cut-off and cv_hold covers that direction,
%   Phi is then held at the cut-off until |I| falls to I_cutoff (a C-rate)
%   or cv_max_time passes. Under constant voltage the currents follow
%   directly from Butler-Volmer at the held potential (REDUCED_RHS).
%
%   Output fields
%     t        time [s]                             (nt x 1)
%     SOC      state of charge                      (nt x 1)
%     Phi      electrode potential [V]              (nt x 1)
%     I        applied (ramped) current [A]; during a CV hold, the
%              resulting current sum_i N_i*A_i*jtr_i (nt x 1)
%     cv       true during constant-voltage holds   (nt x 1)
%     c_av     average concentration [mol/m^3]      (nt x n_groups)
%     theta    c_av/c_m                             (nt x n_groups)
%     jtr      reaction current density [A/m^2]     (nt x n_groups)
%     sig_ext  external radial stress, V_c*sigma/(R_g*T) (nt x n_groups)
%     mu       chemo-mechanical potential [J/mol]   (nt x n_groups)
%     s        phase state (0, +1, -1)              (nt x n_groups)
%     params, protocol, n_events, runtime
%   protocol.t_split holds the realised segment boundaries and
%   protocol.steps one entry per CC or CV step (segment, mode, t_start, t_end).

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
check_initial_state(params);
n        = params.n_groups;

check_stress_regime(params);

if params.verbose
    fprintf('MechoCPM: %d size groups, Lambda = %g (A_int = %.4g), %s\n', ...
        n, params.Lambda, params.A_int, medium_label(params));
    if params.phase_separation
        fprintf('  binodal points [%.4f %.4f], spinodal points [%.4f %.4f]\n', ...
            params.c_b, params.c_sp);
    else
        fprintf('  no phase separation (Lambda >= Lambda_c = %.4g): solid solution\n', ...
            params.Lambda_c);
    end
end

% Initial state: uniform concentration, single-phase particles. The external
% stress is zero, except for an elastic shell formed at a filling c_form
% different from the initial one.
c0      = params.soc_waypoints(1)*params.c_m*ones(n, 1);
s       = zeros(n, 1);
sig_ext = elastic_stress(c0, params);
Phi0    = -mean(mu_hys(c0/params.c_m, s, sig_ext, params))/params.F;
w       = [c0; zeros(n, 1); sig_ext; Phi0];

t_all = []; w_all = []; s_all = []; I_all = []; cv_all = [];
n_events = 0;

% CCCV: protocol events end each constant-current step at the cut-off
% potential, and constant-voltage holds follow (see the header).
cccv    = strcmpi(params.protocol_type, 'CCCV');
weights = params.N_part.*params.Apart;
n_seg   = numel(protocol.I);
t_seg   = protocol.t_split;                  % realised segment boundaries
steps   = struct('segment', {}, 'mode', {}, 't_start', {}, 't_end', {});
I_last  = 0;                                 % current at the end of the last segment
if cccv
    cv_max = params.cv_max_time;
    if isempty(cv_max), cv_max = 3600/params.I_cutoff; end
    I_cut = params.I_cutoff*params.Q/3600;
end

for b = 1:n_seg
    I = protocol.I(b);
    I_prev = I_last;

    if I_prev ~= 0 && sign(I) ~= sign(I_prev)
        params.tau_I = params.tau_I_reversal;
    else
        params.tau_I = params.tau_I_default;
    end

    sig_ext = w(2*n+1:3*n);                  % continuous across blocks
    if cccv
        tspan = t_seg(b) + [0, protocol.duration(b)];
    else
        tspan = protocol.t_split(b:b+1);
    end

    % Constant-current step (ends at the cut-off potential under CCCV)
    step = [];
    if cccv && I ~= 0
        lith = I < 0;
        step = struct('mode', 'CC', 'Phi_cut', params.V_cutoff(2 - lith)/1000, ...
                      'sense', 2*lith - 1, 'max_step', step_size(params, tspan));
    end
    hit = false;
    if cccv && I ~= 0 && step.sense*(w(end) - step.Phi_cut) <= 0
        % already at or beyond the cut-off: no constant-current step
        hit = true;
        t_b = []; w_b = zeros(0, 3*n + 1); s_b = zeros(0, n); ne = 0; I_b = [];
        if params.verbose
            fprintf('  block %d/%d: Phi already at the cut-off %+.1f mV, no CC step\n', ...
                b, n_seg, 1000*step.Phi_cut);
        end
    elseif strcmpi(params.solver, 'ode15s')
        % Reduced ODE in (c_avg, sig_ext): no consistent initial conditions needed
        if params.verbose
            fprintf('  block %d/%d: I = %+.3e A, %.4g h (ode15s%s)\n', ...
                b, n_seg, I, diff(tspan)/3600, cc_label(step));
        end
        [t_b, w_b, s_b, ne, ~, hit, s] = solve_block_ode15s(tspan, w(1:n), sig_ext, w(end), ...
            s, I, I_prev, params, step);
        I_b = ramp_current(t_b, tspan(1), I_prev, I, params.tau_I);
    else
        % Fully implicit DAE in (c_avg, jtr, sig_ext, Phi)
        [w0, wdot0] = consistent_initial_conditions(tspan(1), w(1:n), s, ...
            sig_ext, w(end), I, I_prev, tspan(1), params);
        if params.verbose
            res0 = dae_residual(tspan(1), w0, wdot0, s, I, I_prev, tspan(1), params);
            fprintf('  block %d/%d: I = %+.3e A, %.4g h (ode15i), |res(t0)| = %.1e\n', ...
                b, n_seg, I, diff(tspan)/3600, norm(res0));
        end
        [t_b, w_b, s_b, ne] = solve_block(tspan, w0, wdot0, s, I, I_prev, params);
        s = s_b(end, :)';
        I_b = ramp_current(t_b, tspan(1), I_prev, I, params.tau_I);
    end
    n_events = n_events + ne;
    I_last   = I;

    if ~isempty(t_b)
        if params.verbose
            fprintf('           Phi %+.2f -> %+.2f mV (min %+.2f, max %+.2f), %d events%s\n', ...
                1000*w_b(1, end), 1000*w_b(end, end), 1000*min(w_b(:, end)), ...
                1000*max(w_b(:, end)), ne, cut_note(hit, w_b, t_b, params));
        end
        w = w_b(end, :)';
        t_all  = [t_all; t_b];                       %#ok<AGROW>
        w_all  = [w_all; w_b];                       %#ok<AGROW>
        s_all  = [s_all; s_b];                       %#ok<AGROW>
        I_all  = [I_all; I_b];                       %#ok<AGROW>
        cv_all = [cv_all; false(numel(t_b), 1)];     %#ok<AGROW>
        steps(end+1) = struct('segment', b, 'mode', 'CC', 't_start', t_b(1), 't_end', t_b(end)); %#ok<AGROW>
    end

    % Constant-voltage hold at the cut-off
    if cccv && I ~= 0 && hit && holds(params.cv_hold, lith)
        if isempty(t_all), t0 = tspan(1); else, t0 = t_all(end); end
        tspan_cv = t0 + [0, cv_max];
        cv = struct('mode', 'CV', 'Phi_hold', step.Phi_cut, 'I_cut', I_cut, ...
                    'max_step', step.max_step);
        if params.verbose
            fprintf('           CV hold at %+.1f mV until |I| < C/%g\n', ...
                1000*cv.Phi_hold, 1/params.I_cutoff);
        end
        [t_c, w_c, s_c, ne, ~, done, s] = solve_block_ode15s(tspan_cv, w(1:n), ...
            w(2*n+1:3*n), cv.Phi_hold, s, I, I_prev, params, cv);
        I_c = w_c(:, n+1:2*n)*weights;
        n_events = n_events + ne;
        if ~done
            warning('MechoCPM:cvMaxTime', ...
                ['CV hold of block %d reached cv_max_time = %g s with |I| = C/%.3g ' ...
                 '(I_cutoff = C/%g).'], b, cv_max, params.Q/(3600*abs(I_c(end))), 1/params.I_cutoff);
        end
        if params.verbose
            fprintf('           %.3g h, SoC %+.4f, |I| -> C/%.3g, %d events\n', ...
                (t_c(end) - t_c(1))/3600, soc_change(w_c, params), ...
                params.Q/(3600*abs(I_c(end))), ne);
        end
        w = w_c(end, :)';
        I_last = I_c(end);
        t_all  = [t_all; t_c];                       %#ok<AGROW>
        w_all  = [w_all; w_c];                       %#ok<AGROW>
        s_all  = [s_all; s_c];                       %#ok<AGROW>
        I_all  = [I_all; I_c];                       %#ok<AGROW>
        cv_all = [cv_all; true(numel(t_c), 1)];      %#ok<AGROW>
        steps(end+1) = struct('segment', b, 'mode', 'CV', 't_start', t_c(1), 't_end', t_c(end)); %#ok<AGROW>
    end
    if cccv && ~isempty(t_all)
        t_seg(b + 1) = t_all(end);
    end
end

protocol.t_split = t_seg;
protocol.steps   = steps;
sol = postprocess(t_all, w_all, s_all, I_all, params, protocol);
sol.cv = cv_all;
sol.n_events = n_events;
sol.runtime  = toc(timer);

if params.verbose
    fprintf('  done: %d phase-change events, Phi in [%+.2f, %+.2f] mV, %.1f s\n', ...
        n_events, 1000*min(sol.Phi), 1000*max(sol.Phi), sol.runtime);
end
end

function ok = holds(cv_hold, lith)
% True if a CV hold follows a lithiation (lith) or delithiation step.
switch lower(cv_hold)
    case 'both', ok = true;
    case 'lith', ok = lith;
    case 'delith', ok = ~lith;
    otherwise, ok = false;
end
end

function h = step_size(params, tspan)
% Maximum ode15s step: as in SOLVE_BLOCK_ODE15S, 1/1000 of the segment.
if isempty(params.max_step), h = diff(tspan)/1000; else, h = params.max_step; end
end

function s = cc_label(step)
if isempty(step), s = ''; else, s = sprintf(', CC to %+.1f mV', 1000*step.Phi_cut); end
end

function s = cut_note(hit, w_b, ~, params)
if hit
    n = params.n_groups;
    s = sprintf(', cut-off reached at SoC %.4f', soc_of(w_b(end, 1:n), params));
else
    s = '';
end
end

function x = soc_of(c_row, params)
x = c_row*(params.N_part.*params.F.*params.Vol)/params.Q;
end

function d = soc_change(w, params)
n = params.n_groups;
d = soc_of(w(end, 1:n), params) - soc_of(w(1, 1:n), params);
end

function check_initial_state(params)
% The switching rule (paper, Eq. state_var) assumes every particle starts
% single phase outside the spinodal region. A single-phase particle with
% c_s1 < c_avg < c_s2 is unstable, and the first spinodal point it crosses
% would assign the wrong morphology (e.g. lithiating from 0.5 crosses c_s2
% and would become a lithium-poor shell).
if ~params.phase_separation, return; end
c0 = params.soc_waypoints(1);
if c0 > params.c_sp(1) && c0 < params.c_sp(2)
    error('MechoCPM:unstableStart', ...
        ['The initial state c_avg = %g lies between the spinodal points [%.4f %.4f] ' ...
         '(Lambda = %g), where a single-phase particle is unstable and the ' ...
         'phase-switching rule is not defined. Start outside this range, e.g. ' ...
         'soc_waypoints(1) <= %.4f or >= %.4f.'], ...
        c0, params.c_sp(1), params.c_sp(2), params.Lambda, params.c_sp(1), params.c_sp(2));
end
end

function check_stress_regime(params)
% The plastic-shell model behaves as a rigid-perfectly plastic shell when
% u0 << swelling rate and Lambda*(1+3*nu_bar)*u0 << k*sigma_y. Warn otherwise.
if ~strcmp(external_stress_law(params), 'plastic'), return; end
rate = min(params.c_rate)/3600;                    % filling-fraction rate [1/s]
K = params.Lambda*(1 + 3*params.nu_bar);
if params.u0_used > 0.1*rate
    warning('MechoCPM:u0Regime', ...
        ['u0 = %g 1/s is not small compared with the swelling rate C-rate/3600 = %g 1/s: ' ...
         'the external stress will not reach +-sigma_y. Reduce u0.'], params.u0_used, rate);
end
if K*params.u0_used > 0.1*params.k*params.sigma_y
    warning('MechoCPM:u0Regime', ...
        ['Lambda*(1+3*nu_bar)*u0 = %g is not small compared with k*sigma_y = %g: ' ...
         'the shell will not behave as rigid-plastic. Reduce u0 or increase k.'], ...
        K*params.u0_used, params.k*params.sigma_y);
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
