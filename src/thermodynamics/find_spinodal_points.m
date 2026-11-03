function c_sp = find_spinodal_points(params)
%FIND_SPINODAL_POINTS Spinodal points of a particle with mechanical coupling.
%
%   c_sp = FIND_SPINODAL_POINTS(params) returns the spinodal points
%   c_sp = [c_s1, c_s2] (dimensionless, c_s1 < c_s2): the roots of the
%   coherent-spinodal stability condition (paper, "Stability analysis")
%
%       mu_chem'(c) + 2*A_int = 0,   A_int = Lambda*(1 + 3*nu_bar)/(3*(1 + nu_bar)),
%
%   between which a single-phase particle is unstable. The roots are found
%   by bracketing either side of the global minimum of mu_chem'(c)
%   (MIN_MU_CHEM_PRIME), so the result does not depend on an initial guess.
%   This assumes mu_chem'(c) + 2*A_int changes sign exactly once on each
%   side of that minimum, i.e. a single miscibility gap.
%
%   If 2*A_int > -min mu_chem' (Lambda >= Lambda_c = 1.8 for nu_bar = 0.5
%   and the default mu_chem), mechanics suppresses phase separation: a
%   warning is issued and NaN is returned for both points.
%
%   See also FIND_BINODAL_POINTS, MU_CHEM_PRIME.

K = 2*params.A_int;
stability = @(c) mu_chem_prime(c) + K;

% Most unstable concentration: global minimum of mu_chem'(c) on (0, 1).
edge = 1e-9;
[c_min, dmu_min] = min_mu_chem_prime();
s_min = dmu_min + K;

if s_min >= 0
    warning('find_spinodal_points:stable', ...
        'No spinodal points: mechanical coupling suppresses phase separation (Lambda = %.4g, 2*A_int = %.4g).', params.Lambda, K);
    c_sp = [NaN, NaN];
    return
end

% One root on each side of the minimum; stability(c) > 0 near c = 0 and c = 1.
c_s1 = fzero(stability, [edge, c_min]);
c_s2 = fzero(stability, [c_min, 1 - edge]);
c_sp = [c_s1, c_s2];
end
