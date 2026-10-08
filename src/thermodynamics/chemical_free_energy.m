function [f, mu, dmu] = chemical_free_energy(c)
%CHEMICAL_FREE_ENERGY Homogeneous chemical free energy and its derivatives.
%
%   [f, mu, dmu] = CHEMICAL_FREE_ENERGY(c) returns, elementwise in the
%   normalised concentration c in (0,1), the dimensionless free energy
%   g_chem = f, the chemical potential mu_chem = df/dc (Eq. 72) and its
%   derivative mu_chem' = d^2f/dc^2 (all in units of R_g*T).
%
%   This is the regular-solution free energy with interaction parameter
%   Omega = 3 used in the paper (an LFP-like material):
%
%       f      = c*log(c) + (1-c)*log(1-c) + 3*c*(1-c)
%       mu     = log(c/(1-c)) + 3*(1-2*c)
%       dmu    = 1/(c*(1-c)) - 6
%
%   To use another material, replace the three expressions below. f, mu and
%   dmu are kept in this one file so they stay mutually consistent;
%   tests/test_thermodynamics.m checks that mu = f' and dmu = mu'.

Omega = 3;   % regular-solution interaction parameter [R_g*T]

f   = c.*log(c) + (1 - c).*log(1 - c) + Omega*c.*(1 - c);
mu  = log(c./(1 - c)) + Omega*(1 - 2*c);
dmu = 1./(c.*(1 - c)) - 2*Omega;
end
