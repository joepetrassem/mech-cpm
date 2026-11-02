function [t, w, s_hist, n_events, Phi_last] = solve_block_ode15s(tspan, c0, sig0, Phi0, s, I, I_prev, params)
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

n        = params.n_groups;
t_block  = tspan(1);
t_start  = tspan(1);
y0       = [c0(:); sig0(:)];
Phi_last = Phi0;
t = []; w = []; s_hist = [];
n_events = 0;

while true
    s_now = s;
    opts = odeset('RelTol', params.RelTol, 'AbsTol', params.AbsTol);
    if params.phase_separation
        opts = odeset(opts, 'Events', @(tt, yy) phase_change_events(tt, yy, [], s_now, params));
    end
    lastwarn('');
    if params.phase_separation
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
    s = update_phase_state(s, hits);
    n_events = n_events + numel(hits);
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
        [ydot, Phi_last] = reduced_rhs(tt, yy, s_now, Phi_last, I, I_prev, t_block, params);
    end

    function W = full_state(tt, Y, s_seg)
        % Rebuild [c_av, jtr, sig_ext, Phi] at the output times.
        W = zeros(numel(tt), 3*n + 1);
        for k = 1:numel(tt)
            c_k = Y(k, 1:n)';  sig_k = Y(k, n+1:2*n)';
            [Phi_last, jtr_k] = solve_potential(tt(k), c_k, s_seg, sig_k, Phi_last, ...
                I, I_prev, t_block, params);
            W(k, :) = [c_k; jtr_k; sig_k; Phi_last]';
        end
    end
end
