function [dFdy, dFdyp] = dae_jacobian(t, w, wdot, s, I, I_prev, t_block, params) %#ok<INUSL>
%DAE_JACOBIAN Jacobian of DAE_RESIDUAL for ODE15I.
%
%   [dFdy, dFdyp] = DAE_JACOBIAN(t, w, wdot, s, I, I_prev, t_block, params)
%   returns the sparse partial derivatives of the residual F(t, w, w') with
%   respect to w and w'. State ordering w = [c_av (n); jtr (n); sig_ext (n); Phi].
%
%   - Butler-Volmer rows: derivatives with respect to c_av, sig_ext and Phi
%     by central differences, one particle at a time but vectorised (each
%     particle's current depends only on its own c_av and sig_ext), so any
%     change to MU_HYS or EXCHANGE_CURRENT is picked up automatically.
%   - External-stress rows: analytic, using the derivative returned by
%     SIG_IDEAL. The target stress switches steeply with the surface
%     velocity, and finite differences there mislead Newton's method.
%   - Lithium and current conservation rows: exact.

n   = params.n_groups;
c   = w(1:n);
sig = w(2*n+1:3*n);
Phi = w(3*n+1);
c_dot   = wdot(1:n);
sig_dot = wdot(2*n+1:3*n);
k = params.k;

ic = 1:n;  ij = n+1:2*n;  is = 2*n+1:3*n;  ip = 3*n+1;
N = 3*n + 1;

%% Butler-Volmer derivatives (diagonal in particles)
hc = 1e-6*params.c_m;
hs = 1e-6;
hp = 1e-7;
bv = @(cc, ss, pp) exchange_current(cc, s, pp, ss, params);
dBV_dc   = (bv(c + hc, sig, Phi) - bv(c - hc, sig, Phi))/(2*hc);
dBV_dsig = (bv(c, sig + hs, Phi) - bv(c, sig - hs, Phi))/(2*hs);
dBV_dPhi = (bv(c, sig, Phi + hp) - bv(c, sig, Phi - hp))/(2*hp);

%% External-stress derivatives
switch external_stress_law(params)
    case 'plastic'
        K = params.Lambda*(1 + 3*params.nu_bar);
        [~, dsig_du] = sig_ideal(surface_velocity(sig_dot, c_dot, params), params);
        dF3_dsigdot = 1 - k*dsig_du/K;        % F3 = sig_dot + k*(sig - sig_ideal(u_dot))
        dF3_dcdot   =   - k*dsig_du/params.c_m;
        dF3_dsig    = k*ones(n, 1);
    case 'elastic'                            % F3 = sig_dot + K_eff*c_dot/c_m
        dF3_dsigdot = ones(n, 1);
        dF3_dcdot   = params.K_eff/params.c_m*ones(n, 1);
        dF3_dsig    = zeros(n, 1);
    otherwise                                 % F3 = sig_dot + k*sig
        dF3_dsigdot = ones(n, 1);
        dF3_dcdot   = zeros(n, 1);
        dF3_dsig    = k*ones(n, 1);
end

%% Assemble
rows = [ic, ij, ij, ij, ij, is, ip*ones(1, n)];
cols = [ij, ic, ij, is, ip*ones(1, n), is, ij];
vals = [(params.Apart./(params.F*params.Vol))', dBV_dc', -ones(1, n), dBV_dsig', dBV_dPhi', ...
        dF3_dsig', -(params.N_part.*params.Apart)'];
dFdy = sparse(rows, cols, vals, N, N);

rows = [ic, is, is];
cols = [ic, ic, is];
vals = [ones(1, n), dF3_dcdot', dF3_dsigdot'];
dFdyp = sparse(rows, cols, vals, N, N);
end
