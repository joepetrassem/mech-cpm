function [ydot, Phi, jtr] = reduced_rhs(t, y, s, Phi_guess, I, I_prev, t_block, params)
%REDUCED_RHS Right-hand side of the reduced ODE in (c_avg, sig_ext).
%
%   [ydot, Phi, jtr] = REDUCED_RHS(t, y, s, Phi_guess, I, I_prev, t_block, params)
%   with y = [c_av (n); sig_ext (n)] returns ydot = [dc_av/dt; dsig_ext/dt].
%
%   The algebraic variables and the implicit stress law of the DAE are
%   solved inside this function, so that the integrator (ODE15S) sees an
%   ordinary ODE:
%     1. Phi from current conservation, sum_i N_i*A_i*jtr_i(Phi) = I_ramp(t)
%        (SOLVE_POTENTIAL; unique, since the total current increases with Phi);
%     2. jtr_i = BV(c_i, s_i, Phi, sig_i), and dc_i/dt = -jtr_i*A_i/(F*V_i);
%     3. dsig_i/dt from the implicit law dsig/dt = -k*(sig - sig_ideal(u_dot))
%        (SOLVE_STRESS_RATE; unique root, found by bisection).
%   These are exactly the equations of DAE_RESIDUAL, so the two formulations
%   describe the same model.
%
%   Phi_guess warm-starts the root finder (e.g. the Phi of the previous call).
%
%   See also SOLVE_BLOCK_ODE15S, DAE_RESIDUAL.

n       = params.n_groups;
c_av    = y(1:n);
sig_ext = y(n+1:2*n);

[Phi, jtr] = solve_potential(t, c_av, s, sig_ext, Phi_guess, I, I_prev, t_block, params);
c_dot   = -jtr.*params.Apart./(params.F*params.Vol);
sig_dot = solve_stress_rate(sig_ext, c_dot, params);
ydot    = [c_dot; sig_dot];
end
