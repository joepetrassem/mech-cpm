function [c_min, dmu_min] = min_mu_chem_prime()
%MIN_MU_CHEM_PRIME Global minimum of mu_chem'(c) on (0, 1).
%
%   [c_min, dmu_min] = MIN_MU_CHEM_PRIME() returns the concentration c_min
%   of the most unstable state of the homogeneous chemical potential and the
%   value dmu_min = mu_chem'(c_min). The minimum is located on a grid of
%   2001 points (so that a mu_chem' with several local minima is handled),
%   then refined with FMINBND between the neighbouring grid points.
%
%   Phase separation is possible only if 2*A_int < -dmu_min, which defines
%   the critical coupling Lambda_c (FINALISE_PARAMETERS).
c  = linspace(1e-4, 1 - 1e-4, 2001);
d  = mu_chem_prime(c);
[~, k] = min(d);
lo = c(max(k - 1, 1));  hi = c(min(k + 1, numel(c)));
[c_min, dmu_min] = fminbnd(@mu_chem_prime, lo, hi, optimset('TolX', 1e-12));
end
