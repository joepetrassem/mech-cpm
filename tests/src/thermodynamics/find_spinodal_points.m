function c_sp = find_spinodal_points(params)
%FIND_SPINODAL_POINTS Spinodal points of a particle with mechanical coupling.
%
%   c_sp = FIND_SPINODAL_POINTS(params) returns the spinodal points
%   c_sp = [c_s1, c_s2] (dimensionless, c_s1 < c_s2): the roots of the
%   stability condition (Eq. 71)
%
%       mu_chem'(c) + K = 0,   K = 2*gamma*alpha*(1 + 3*nu_bar)/(3*(1 + nu_bar)),
%
%   between which a single-phase particle is unstable. The roots are found
%   by bracketing either side of the minimum of mu_chem'(c), so the result
%   does not depend on an initial guess.
%
%   If K is large enough that mu_chem'(c) + K > 0 everywhere, mechanics
%   suppresses phase separation: a warning is issued and NaN is returned
%   for both points.
%
%   See also FIND_BINODAL_POINTS, MU_CHEM_PRIME.

cons = 1 + 3*params.nu_bar;                          % 1 + 3*nu_bar
K = 2*params.gamma*params.alpha*cons/(cons + 2);     % cons + 2 = 3*(1 + nu_bar)
stability = @(c) mu_chem_prime(c) + K;

% Most unstable concentration: minimum of mu_chem'(c) on (0, 1).
edge = 1e-9;
[c_min, s_min] = fminbnd(stability, edge, 1 - edge, optimset('TolX', 1e-12));

if s_min >= 0
    warning('find_spinodal_points:stable', ...
        'No spinodal points: mechanical coupling suppresses phase separation (K = %.4g).', K);
    c_sp = [NaN, NaN];
    return
end

% One root on each side of the minimum; stability(c) > 0 near c = 0 and c = 1.
c_s1 = fzero(stability, [edge, c_min]);
c_s2 = fzero(stability, [c_min, 1 - edge]);
c_sp = [c_s1, c_s2];
end
