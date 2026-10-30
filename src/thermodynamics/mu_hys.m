function mu = mu_hys(c_avg, s, sig_ext, params)
%MU_HYS Surface chemo-mechanical potential of a particle [J/mol].
%
%   mu = MU_HYS(c_avg, s, sig_ext, params) evaluates, elementwise, the
%   chemo-mechanical potential at the particle surface for a particle with
%   dimensionless average concentration c_avg (= c/c_m), phase state s and
%   external radial stress sig_ext. s and sig_ext may be scalars or arrays
%   the same size as c_avg. The result is R_g*T times the dimensionless
%   potential, so that eta = (F*Phi + mu)/(R_g*T).
%
%   Phase states
%     s =  0 : single phase, c_c = c_s = c_avg
%     s = +1 : core-shell, lithium-rich shell:  c_c = c_b1, c_s = c_b2
%     s = -1 : core-shell, lithium-poor shell:  c_c = c_b2, c_s = c_b1
%
%   With stresses scaled by R_g*T/V_c (paper, "Non-dimensionalisation"),
%   the surface potential is (paper: mu_single and mu_sep)
%       single phase : mu = mu_chem(c_avg) - sig_ext
%       core-shell   : mu = mu_chem(c_s) - 2*A_int*(c_avg - c_s) - sig_ext
%   with A_int = Lambda*(1 + 3*nu_bar)/(3*(1 + nu_bar)). The external stress
%   translates both branches rigidly; Lambda (through A_int) tilts the
%   mixed-phase branch.
%
%   See also EXCHANGE_CURRENT, FIND_BINODAL_POINTS.

sz      = size(c_avg);
s       = expand(s, sz);
sig_ext = expand(sig_ext, sz);

mu = mu_single_phase(c_avg, sig_ext, params);

c_b = sort(params.c_b);
for state = [1, -1]
    m = (s == state);
    if ~any(m(:)), continue; end
    if state == 1
        c_s = c_b(2);  c_c = c_b(1);     % lithium-rich shell
    else
        c_s = c_b(1);  c_c = c_b(2);     % lithium-poor shell
    end
    mu(m) = mu_core_shell(c_avg(m), sig_ext(m), c_c, c_s, params);
end

mu = params.R_gT*mu;
end

% -------------------------------------------------------------------------
function mu = mu_single_phase(c_avg, sig_ext, params)
% Dimensionless potential of a single-phase particle (sigma_kk = 3*sig_ext).
mu = mu_chem(c_avg) - sig_ext;
end

function mu = mu_core_shell(c_avg, sig_ext, ~, c_s, params)
% Dimensionless surface (shell) potential of a core-shell particle.
% c_c is not needed explicitly: c_avg - c_s = psi*(c_c - c_s) by the lever rule.
mu = mu_chem(c_s) - 2*params.A_int*(c_avg - c_s) - sig_ext;
end

function x = expand(x, sz)
if isscalar(x)
    x = x*ones(sz);
else
    x = reshape(x, sz);
end
end
