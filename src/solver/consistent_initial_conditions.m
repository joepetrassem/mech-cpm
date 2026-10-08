function [w0, wdot0] = consistent_initial_conditions(t0, c_av, s, sig_ext, Phi_guess, I, I_prev, t_block, params)
%CONSISTENT_INITIAL_CONDITIONS Consistent (w, w') for ODE15I.
%
%   Given the differential states c_av and sig_ext and the phase states s,
%   solves the algebraic constraints (Butler-Volmer and conservation of
%   current) for Phi and jtr at time t0, then evaluates the time
%   derivatives. Used at t = 0, at every block start and after every
%   phase-change event.
%
%   Phi is found by bracketing + FZERO on the scalar total-current
%   equation, which is monotonically increasing in Phi.

c_av = c_av(:);  s = s(:);  sig_ext = sig_ext(:);

[Phi, jtr] = solve_algebraic(t0, c_av, s, sig_ext, Phi_guess, I, I_prev, t_block, params);

c_dot   = -jtr.*params.Apart./(params.F*params.Vol);
sig_dot = -params.k*(sig_ext - sig_ideal(I, params));

% Derivatives of the algebraic variables by a forward difference. They are
% not needed for consistency but give ODE15I a better first step.
dt = 1e-3;
[Phi2, jtr2] = solve_algebraic(t0 + dt, c_av + dt*c_dot, s, sig_ext + dt*sig_dot, ...
                               Phi, I, I_prev, t_block, params);

w0    = [c_av; jtr; sig_ext; Phi];
wdot0 = [c_dot; (jtr2 - jtr)/dt; sig_dot; (Phi2 - Phi)/dt];
end

function [Phi, jtr] = solve_algebraic(t, c_av, s, sig_ext, Phi_guess, I, I_prev, t_block, params)
I_now   = ramp_current(t, t_block, I_prev, I, params.tau_I);
weights = params.N_part.*params.Apart;
g = @(Phi) sum(weights.*exchange_current(c_av, s, Phi, sig_ext, params)) - I_now;
Phi = bracket_and_solve(g, Phi_guess);
jtr = exchange_current(c_av, s, Phi, sig_ext, params);
end

function x = bracket_and_solve(g, x0)
g0 = g(x0);
if g0 == 0
    x = x0;
    return
end
direction = -sign(g0);        % g is increasing in Phi
step = 1e-3;                  % [V]
x1 = x0; g1 = g0;
for k = 1:60
    x1 = x0 + direction*step;
    g1 = g(x1);
    if ~isfinite(g1)
        error('MechoCPM:bracket', 'Non-finite current while bracketing Phi.');
    end
    if sign(g1) ~= sign(g0), break; end
    x0 = x1; g0 = g1;
    step = 2*step;
end
if sign(g1) == sign(g0)
    error('MechoCPM:bracket', 'Could not bracket the electrode potential.');
end
x = fzero(g, sort([x0, x1]), optimset('TolX', 1e-14));
end
