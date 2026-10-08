function [t, w, s_hist, n_events] = solve_block(tspan, w0, wdot0, s, I, I_prev, params)
%SOLVE_BLOCK Integrate one constant-current block with ODE15I.
%
%   Integrates until either the end of the block or a phase-change event.
%   At each event the phase state of the affected particle(s) is updated,
%   consistent initial conditions are recomputed and integration resumes.
%
%   Outputs: time vector t, states w (one row per time), phase states
%   s_hist (one row per time) and the number of events handled.

n        = params.n_groups;
t_block  = tspan(1);
t_start  = tspan(1);
t = []; w = []; s_hist = [];
n_events = 0;

while true
    s_now = s;   % frozen copy captured by the function handles below
    opts = odeset('RelTol', params.RelTol, 'AbsTol', params.AbsTol, ...
        'Events', @(tt, ww, yp) phase_change_events(tt, ww, yp, s_now, params));
    [t_p, w_p, te, ~, ie] = ode15i( ...
        @(tt, ww, yp) dae_residual(tt, ww, yp, s_now, I, I_prev, t_block, params), ...
        [t_start, tspan(2)], w0, wdot0, opts);

    t      = [t; t_p];                                   %#ok<AGROW>
    w      = [w; w_p];                                   %#ok<AGROW>
    s_hist = [s_hist; repmat(s_now', numel(t_p), 1)];    %#ok<AGROW>

    if t_p(end) >= tspan(2)
        break
    end
    if isempty(ie)
        error('MechoCPM:solverStopped', ...
            'ode15i stopped at t = %.6g s (block ends at %.6g s) without an event.', ...
            t_p(end), tspan(2));
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
    w_end   = w_p(end, :)';
    [w0, wdot0] = consistent_initial_conditions(t_start, w_end(1:n), s, ...
        w_end(2*n+1:3*n), w_end(end), I, I_prev, t_block, params);
end
end

function s = update_phase_state(s, event_ids)
% homogeneous -> core-shell: crossing spin point 1 gives s = +1 (lithiating
% morphology), spin point 2 gives s = -1 (delithiating morphology).
% core-shell -> homogeneous: reaching either star point.
for e = event_ids(:)'
    p     = ceil(e/2);
    which = e - 2*(p - 1);
    if s(p) == 0
        if which == 1
            s(p) = 1;
        else
            s(p) = -1;
        end
    else
        s(p) = 0;
    end
end
end
