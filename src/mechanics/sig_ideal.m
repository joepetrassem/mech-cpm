function sig = sig_ideal(u_dot, params)
%SIG_IDEAL Target external stress of a rapidly yielding plastic shell.
%
%   sig = SIG_IDEAL(u_dot, params) = -sigma_y*tanh(u_dot/u0)
%
%   u_dot is the dimensionless surface velocity of each particle [1/s]
%   (see SURFACE_VELOCITY). The shell opposes the motion of the surface: an
%   expanding particle (u_dot > 0, lithiation) is compressed towards
%   -sigma_y, a contracting one is pulled outwards towards +sigma_y.
%   sigma_y is the scaled yield stress V_c*sigma_y/(R_g*T); u0 is a
%   reference surface velocity [1/s] that sets how sharply the target
%   switches between -+sigma_y.
%
%   To model another external medium, change this function, SURFACE_VELOCITY
%   and the relaxation law in DAE_RESIDUAL (and SOLVE_STRESS_RATE, which
%   solves that law for the stress rate).
%
%   See also SURFACE_VELOCITY, SOLVE_STRESS_RATE, DAE_RESIDUAL.
sig = -params.sigma_y*tanh(u_dot/params.u0);
end
