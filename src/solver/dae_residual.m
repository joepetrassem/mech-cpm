function res = dae_residual(t, w, wdot, s, I, I_prev, t_block, params)
%DAE_RESIDUAL Fully implicit residual F(t, w, w') = 0 for ODE15I.
%
%   State vector  w = [c_av (n); jtr (n); sig_ext (n); Phi (1)]
%
%   Equations, for each size group i:
%     dc_i/dt     = -jtr_i*A_i/(F*V_i)                 (Li conservation, Eq. 74)
%     0           = BV(c_i, s_i, Phi, sig_i) - jtr_i    (Butler-Volmer, Eq. 76)
%     dsig_i/dt   = -k*(sig_i - sig_ideal(I))           (external stress, Eq. 78)
%   and globally
%     0           = I_ramp(t) - sum_i N_i*A_i*jtr_i     (current conservation, Eq. 77)

n = params.n_groups;

c_av    = w(1:n);
jtr     = w(n+1:2*n);
sig_ext = w(2*n+1:3*n);
Phi     = w(3*n+1);

c_dot   = wdot(1:n);
sig_dot = wdot(2*n+1:3*n);

if any(~isfinite(w))
    error('MechoCPM:nonFinite', 'Non-finite state at t = %g s.', t);
end

res = zeros(3*n + 1, 1);
res(1:n)       = c_dot + jtr.*params.Apart./(params.F*params.Vol);
res(n+1:2*n)   = exchange_current(c_av, s, Phi, sig_ext, params) - jtr;
res(2*n+1:3*n) = sig_dot + params.k*(sig_ext - sig_ideal(I, params));
res(end)       = ramp_current(t, t_block, I_prev, I, params.tau_I) ...
                 - sum(params.N_part.*params.Apart.*jtr);
end
