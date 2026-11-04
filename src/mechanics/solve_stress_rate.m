function sig_dot = solve_stress_rate(sig_ext, c_dot, params)
%SOLVE_STRESS_RATE Stress rate that satisfies the external-stress law.
%
%   sig_dot = SOLVE_STRESS_RATE(sig_ext, c_dot, params) solves, for each
%   particle, the implicit relaxation law
%
%       sig_dot = -k*(sig_ext - sig_ideal(u_dot)),
%       u_dot   = sig_dot/(Lambda*(1 + 3*nu_bar)) + c_dot/c_m,
%
%   for sig_dot. The residual  g(sig_dot) = sig_dot + k*(sig_ext +
%   sigma_y*tanh(u_dot/u0))  increases strictly with sig_dot and changes sign
%   on [-k*(sig_ext + sigma_y), -k*(sig_ext - sigma_y)], so the root is
%   unique and is found by vectorised bisection to machine precision.
%   Used for consistent initial conditions; inside the time integration
%   ODE15I solves the same equation as part of the DAE.
%
%   For Lambda = 0 the particle never loads the shell and sig_dot = -k*sig_ext
%   (the stress stays at zero). For the elastic shell (external_medium =
%   'elastic') the rate is explicit: sig_dot = -K_eff*c_dot/c_m.
sig_ext = sig_ext(:);  c_dot = c_dot(:);
k = params.k;
switch external_stress_law(params)
    case 'zero'
        sig_dot = -k*sig_ext;
        return
    case 'elastic'
        sig_dot = -params.K_eff*c_dot/params.c_m;      % explicit for the elastic shell
        return
end
g  = @(sd) sd + k*(sig_ext - sig_ideal(surface_velocity(sd, c_dot, params), params));
lo = -k*(sig_ext + params.sigma_y);
hi = -k*(sig_ext - params.sigma_y);
for it = 1:200
    mid = (lo + hi)/2;
    neg = g(mid) < 0;
    lo(neg)  = mid(neg);
    hi(~neg) = mid(~neg);
    if all(hi - lo <= 4*eps(max(max(abs(lo), abs(hi)), realmin))), break; end
end
sig_dot = (lo + hi)/2;
end
