function [w0, wdot0] = consistent_initial_conditions(t0, c_av, s, sig_ext, Phi_guess, I, I_prev, t_block, params)
%CONSISTENT_INITIAL_CONDITIONS Consistent (w, w') for ODE15I.
%
%   Given the differential states c_av and sig_ext and the phase states s,
%   solves the algebraic constraints (Butler-Volmer and conservation of
%   current) for Phi and jtr at time t0, then evaluates the time
%   derivatives. Used at t = 0, at every block start and after every
%   phase-change event.
%
%   Phi and jtr come from SOLVE_POTENTIAL, dsig_ext/dt from
%   SOLVE_STRESS_RATE.

c_av = c_av(:);  s = s(:);  sig_ext = sig_ext(:);

[Phi, jtr] = solve_potential(t0, c_av, s, sig_ext, Phi_guess, I, I_prev, t_block, params);

c_dot   = -jtr.*params.Apart./(params.F*params.Vol);
sig_dot = solve_stress_rate(sig_ext, c_dot, params);   % implicit stress law

% Derivatives of the algebraic variables by a forward difference. They are
% not needed for consistency but give ODE15I a better first step.
dt = 1e-3;
[Phi2, jtr2] = solve_potential(t0 + dt, c_av + dt*c_dot, s, sig_ext + dt*sig_dot, ...
                               Phi, I, I_prev, t_block, params);

w0    = [c_av; jtr; sig_ext; Phi];
wdot0 = [c_dot; (jtr2 - jtr)/dt; sig_dot; (Phi2 - Phi)/dt];
end
