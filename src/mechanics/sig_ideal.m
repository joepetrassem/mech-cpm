function [sig, dsig_du] = sig_ideal(u_dot, params)
%SIG_IDEAL Target external stress of a rapidly yielding plastic shell.
%
%   sig = SIG_IDEAL(u_dot, params) = -sigma_y*tanh(u_dot/u0)
%   [sig, dsig_du] = SIG_IDEAL(...) also returns d(sig)/d(u_dot), which the
%   DAE Jacobian (DAE_JACOBIAN) uses.
%
%   u_dot is the dimensionless surface velocity of each particle [1/s]
%   (see SURFACE_VELOCITY). The shell opposes the motion of the surface: an
%   expanding particle (u_dot > 0, lithiation) is compressed towards
%   -sigma_y, a contracting one is pulled outwards towards +sigma_y.
%   sigma_y is the scaled yield stress V_c*sigma_y/(R_g*T); u0 is a
%   reference surface velocity [1/s] that sets how sharply the target
%   switches between -+sigma_y. The value used is params.u0_used: u0 if
%   given, otherwise 5% of the nominal swelling rate min(c_rate)/3600.
%
%   To model another external medium, change this function (both outputs),
%   and if needed SURFACE_VELOCITY, the relaxation law in DAE_RESIDUAL /
%   DAE_JACOBIAN and SOLVE_STRESS_RATE.
%
%   See also SURFACE_VELOCITY, SOLVE_STRESS_RATE, DAE_RESIDUAL, DAE_JACOBIAN.
x   = u_dot/params.u0_used;
sig = -params.sigma_y*tanh(x);
if nargout > 1
    dsig_du = -params.sigma_y*sech(x).^2/params.u0_used;
end
end
