function [Phi, jtr] = solve_potential(t, c_av, s, sig_ext, Phi_guess, I, I_prev, t_block, params)
%SOLVE_POTENTIAL Electrode potential and reaction currents from the state.
%
%   [Phi, jtr] = SOLVE_POTENTIAL(t, c_av, s, sig_ext, Phi_guess, I, I_prev, t_block, params)
%   solves the current-conservation constraint
%
%       sum_i N_i*A_i*jtr_i(Phi) = I_ramp(t),   jtr_i = BV(c_i, s_i, Phi, sig_i),
%
%   for the electrode potential Phi [V], and returns the reaction current
%   densities jtr [A/m^2] at that potential. The total current increases
%   strictly with Phi, so the root is unique: it is bracketed by stepping
%   outwards from Phi_guess and refined with FZERO.
%
%   Used by CONSISTENT_INITIAL_CONDITIONS (ode15i) and REDUCED_RHS (ode15s).

c_av = c_av(:);  s = s(:);  sig_ext = sig_ext(:);
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
step = 1e-4;                  % [V]
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
