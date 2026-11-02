function u_dot = surface_velocity(sig_dot, c_dot, params)
%SURFACE_VELOCITY Dimensionless velocity of the particle surface [1/s].
%
%   u_dot = SURFACE_VELOCITY(sig_dot, c_dot, params) with sig_dot the rate
%   of the scaled external stress [1/s] and c_dot the rate of the average
%   concentration [mol/(m^3 s)]:
%
%       u_dot = sig_dot/(Lambda*(1 + 3*nu_bar)) + c_dot/c_m
%
%   The time derivative of u_surf = sig_ext/(Lambda*(1 + 3*nu_bar)) + c_avg - c_ref:
%   the first term is the elastic compliance of the particle, the second
%   its chemical swelling. Requires Lambda > 0.
K = params.Lambda*(1 + 3*params.nu_bar);       % bulk stiffness of the particle
u_dot = sig_dot/K + c_dot/params.c_m;
end
