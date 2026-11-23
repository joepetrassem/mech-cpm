function [t, w, s_hist, n_events, Phi_last, hit_cutoff, s] = solve_block_ode15s(tspan, c0, sig0, Phi0, s, I, I_prev, params, step)
%SOLVE_BLOCK_ODE15S Integrate one constant-current block as a reduced ODE with ODE15S.
%
%   [t, w, s_hist, n_events, Phi_last] = SOLVE_BLOCK_ODE15S(tspan, c0, sig0, Phi0, s, I, I_prev, params)
%   integrates y = [c_av; sig_ext] with REDUCED_RHS from tspan(1) to
%   tspan(2), starting from c0, sig0 and phase states s. Phi0 is a starting
%   guess for the electrode potential.
%
%   Phase changes are located with PHASE_CHANGE_EVENTS; at each event the
%   phase state is updated (UPDATE_PHASE_STATE) and integration resumes from
%   the same state. No consistent initial conditions are needed: Phi and jtr
%   are not integrated but recomputed from (c_av, sig_ext) whenever needed.
%
%   Outputs have the same layout as SOLVE_BLOCK: w has one row per output
%   time, w = [c_av, jtr, sig_ext, Phi], with jtr and Phi reconstructed at
%   the output times; s_hist holds the phase states; n_events counts the
%   events handled; Phi_last is the final potential (warm start for the next
%   block).
%
%   SOLVE_BLOCK_ODE15S(..., params, step) adds a terminal protocol event for
%   CCCV cycling (see RUN_SIMULATION):
%     step.mode = 'CC', step.Phi_cut [V], step.sense (+1 lithiation, -1
%       delithiation): constant current until Phi reaches Phi_cut;
%     step.mode = 'CV', step.Phi_hold [V], step.I_cut [A]: Phi held at
%       Phi_hold until |I| = |sum_i N_i*A_i*jtr_i| falls to I_cut.
%   step.max_step [s], if present, overrides the default maximum step.
%   hit_cutoff is true if the block ended at the protocol event rather than
%   at tspan(2), and s is the phase state at the end of the block.

if nargin < 9 || isempty(step)
    step = struct('mode', 'CC', 'Phi_cut', []);
end
is_cv   = strcmp(step.mode, 'CV');
has_cut = is_cv || ~isempty(step.Phi_cut);
if is_cv, Phi0 = step.Phi_hold; end
weights = params.N_part.*params.Apart;

n        = params.n_groups;
i_cut    = 2*n*params.phase_separation + 1;    % index of the protocol event
hit_cutoff = false;
t_block  = tspan(1);
t_start  = tspan(1);
y0       = [c0(:); sig0(:)];
Phi_last = Phi0;
t = []; w = []; s_hist = [];
n_events = 0;

while true
    s_now = s;
    if isfield(step, 'max_step') && ~isempty(step.max_step)
        max_step = step.max_step;
    elseif isempty(params.max_step)
        max_step = diff(tspan)/1000;           % resolve the voltage curve in the output
    else
        max_step = params.max_step;
    end
    opts = odeset('RelTol', params.RelTol, 'AbsTol', params.AbsTol, 'MaxStep', max_step);
    use_events = params.phase_separation || has_cut;
    if use_events
        opts = odeset(opts, 'Events', @events);
    end
    lastwarn('');
    if use_events
        [t_p, y_p, te, ~, ie] = ode15s(@rhs, [t_start, tspan(2)], y0, opts);
    else
        [t_p, y_p] = ode15s(@rhs, [t_start, tspan(2)], y0, opts);
        te = [];  ie = [];
    end

    w_p = full_state(t_p, y_p, s_now);
    t      = [t; t_p];                                   %#ok<AGROW>
    w      = [w; w_p];                                   %#ok<AGROW>
    s_hist = [s_hist; repmat(s_now', numel(t_p), 1)];    %#ok<AGROW>

    if t_p(end) >= tspan(2)
        break
    end
    if isempty(ie)
        yl  = y_p(end, :)';
        msg = lastwarn;
        if isempty(msg), msg = '(no warning)'; end
        error('MechoCPM:solverStopped', ...
            ['ode15s stopped at t = %.6g s (block ends at %.6g s, started at %.6g s) without an event.\n' ...
             '  ode15s: %s\n  c_avg/c_m in [%.4f, %.4f]; sig_ext in [%.4f, %.4f]; phase states %s'], ...
            t_p(end), tspan(2), t_block, msg, min(yl(1:n))/params.c_m, max(yl(1:n))/params.c_m, ...
            min(yl(n+1:2*n)), max(yl(n+1:2*n)), mat2str(s_now'));
    end

    % Handle every event that fired at the stopping time
    hits = ie(abs(te - t_p(end)) <= 10*eps(max(1, abs(t_p(end)))));
    if isempty(hits), hits = ie(end); end
    cut  = has_cut && any(hits == i_cut);
    hits = hits(hits ~= i_cut);
    if ~isempty(hits)
        s = update_phase_state(s, hits);
    end
    n_events = n_events + numel(hits);
    if cut
        hit_cutoff = true;
        break
    end
    if n_events > params.max_events_per_block
        error('MechoCPM:tooManyEvents', ...
            'More than %d events in one block (event chattering?).', params.max_events_per_block);
    end

    t_start = t_p(end);
    y0      = y_p(end, :)';
end

    % ---------------------------------------------------------------------
    function ydot = rhs(tt, yy)
        % Nested so that the last Phi warm-starts the next root solve.
        [ydot, Phi_last] = reduced_rhs(tt, yy, s_now, Phi_last, I, I_prev, t_block, params, step);
    end

    function [value, isterminal, direction] = events(tt, yy)
        % Phase-change events, then the protocol event (if any).
        if params.phase_separation
            [value, isterminal, direction] = phase_change_events(tt, yy, [], s_now, params);
        else
            value = zeros(0, 1);  isterminal = value;  direction = value;
        end
        if has_cut
            value      = [value; cutoff_value(tt, yy)];
            isterminal = [isterminal; 1];
            direction  = [direction; -1];
        end
    end

    function v = cutoff_value(tt, yy)
        % Positive until the cut-off is reached, then crosses zero downwards.
        c_k = yy(1:n);  sig_k = yy(n+1:2*n);
        if is_cv
            jtr_k = exchange_current(c_k, s_now, step.Phi_hold, sig_k, params);
            v = abs(sum(weights.*jtr_k)) - step.I_cut;
        else
            Phi_k = solve_potential(tt, c_k, s_now, sig_k, Phi_last, I, I_prev, t_block, params);
            v = step.sense*(Phi_k - step.Phi_cut);
        end
    end

    function W = full_state(tt, Y, s_seg)
        % Rebuild [c_av, jtr, sig_ext, Phi] at the output times.
        W = zeros(numel(tt), 3*n + 1);
        for k = 1:numel(tt)
            c_k = Y(k, 1:n)';  sig_k = Y(k, n+1:2*n)';
            if is_cv
                Phi_last = step.Phi_hold;
                jtr_k = exchange_current(c_k, s_seg, Phi_last, sig_k, params);
            else
                [Phi_last, jtr_k] = solve_potential(tt(k), c_k, s_seg, sig_k, Phi_last, ...
                    I, I_prev, t_block, params);
            end
            W(k, :) = [c_k; jtr_k; sig_k; Phi_last]';
        end
    end
end
