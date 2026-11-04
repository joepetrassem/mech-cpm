function res = dae_residual(t, w, wdot, s, I, I_prev, t_block, params)
%DAE_RESIDUAL Fully implicit residual F(t, w, w') = 0 for ODE15I.
%
%   State vector  w = [c_av (n); jtr (n); sig_ext (n); Phi (1)]
%
%   Equations, for each size group i:
%     dc_i/dt     = -jtr_i*A_i/(F*V_i)                 (Li conservation)
%     0           = BV(c_i, s_i, Phi, sig_i) - jtr_i    (Butler-Volmer)
%     dsig_i/dt   = -k*(sig_i - sig_ideal(u_dot_i))     (external stress)
%       with the surface velocity
%       u_dot_i     = (dsig_i/dt)/(Lambda*(1 + 3*nu_bar)) + (dc_i/dt)/c_m,
%       which makes this equation implicit in dsig_i/dt. For Lambda = 0 (or
%       sigma_y = 0) the particle never loads the shell: dsig_i/dt = -k*sig_i,
%       so the stress stays at zero. For the elastic shell the row is
%       dsig_i/dt = -K_eff*(dc_i/dt)/c_m (see ELASTIC_STRESS).
%   and globally
%     0           = I_ramp(t) - sum_i N_i*A_i*jtr_i     (current conservation)

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
switch external_stress_law(params)
    case 'plastic'
        u_dot = surface_velocity(sig_dot, c_dot, params);
        res(2*n+1:3*n) = sig_dot + params.k*(sig_ext - sig_ideal(u_dot, params));
    case 'elastic'
        res(2*n+1:3*n) = sig_dot + params.K_eff*c_dot/params.c_m;
    otherwise
        res(2*n+1:3*n) = sig_dot + params.k*sig_ext;
end
res(end)       = ramp_current(t, t_block, I_prev, I, params.tau_I) ...
                 - sum(params.N_part.*params.Apart.*jtr);
end
