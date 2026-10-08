function jtr = exchange_current(c_av, s, Phi, sig_ext, params)
%EXCHANGE_CURRENT Butler-Volmer reaction current density per particle [A/m^2].
%
%   jtr = EXCHANGE_CURRENT(c_av, s, Phi, sig_ext, params) with c_av the
%   average concentrations [mol/m^3], s the phase states, Phi the electrode
%   potential [V] and sig_ext the external radial stresses. Returns a
%   column vector.
%
%       eta = (F*Phi + mu)/(R_g*T)                        (Eq. 75)
%       j   = 2*I0*sqrt(theta_surf*(1 - theta_surf))*sinh(eta/2)   (Eq. 76)
%
%   theta_surf = c_surf/c_m is the surface concentration: c_avg/c_m for a
%   single-phase particle (s = 0), and the shell binodal point for a
%   core-shell one (c_b2 for s = +1, c_b1 for s = -1).
%
%   See also MU_HYS.

c_b = sort(params.c_b);
theta = c_av(:)/params.c_m;
s = s(:);

mu  = mu_hys(theta, s, sig_ext, params);        % [J/mol]
eta = (params.F*Phi + mu)/params.R_gT;          % dimensionless overpotential

theta_surf = theta;
theta_surf(s ==  1) = c_b(2);
theta_surf(s == -1) = c_b(1);

jtr = 2*params.I0*sqrt(theta_surf.*(1 - theta_surf)).*sinh(eta/2);
jtr = real(jtr(:));
end
