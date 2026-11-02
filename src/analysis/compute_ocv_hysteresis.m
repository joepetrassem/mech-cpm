function ocv = compute_ocv_hysteresis(params, sig_ext, n_points)
%COMPUTE_OCV_HYSTERESIS Quasi-static single-particle potential (paper Fig. 3).
%
%   ocv = COMPUTE_OCV_HYSTERESIS(params, sig_ext, n_points) evaluates MU_HYS
%   for one particle cycled quasi-statically, without time integration, for
%   one or more external stresses sig_ext.
%
%   Solution branches (paper Fig. 3a,b)
%     single phase : mu_single = mu_chem(c_avg) - sig_ext, for all c_avg.
%                    Stable outside [c_b1, c_b2], metastable on
%                    [c_b1, c_s1] and [c_s2, c_b2], unstable on [c_s1, c_s2].
%     mixed phase  : mu_sep = mu_chem(c_s) - 2*A_int*(c_avg - c_s) - sig_ext on
%                    [c_b1, c_b2]. The rich-shell (s = +1) and poor-shell
%                    (s = -1) forms give the same line, because the core and
%                    shell potentials are equal at equilibrium.
%
%   Cycling paths (paper Fig. 3c,d), in the direction of travel
%     lithiation   : single phase from c_lo to c_s1 (metastable beyond c_b1),
%                    jump onto the mixed-phase line at c_s1 (s: 0 -> +1),
%                    mixed phase to c_b2, single phase to c_hi.
%     delithiation : single phase from c_hi to c_s2 (metastable below c_b2),
%                    jump onto the mixed-phase line at c_s2 (s: 0 -> -1),
%                    mixed phase to c_b1, single phase to c_lo.
%   The jump points appear twice in the path (before and after the switch).
%
%   Without phase separation (Lambda >= Lambda_c, c_b = c_sp = NaN) every
%   state is single phase and stable, there is no mixed-phase line, and
%   both paths follow mu_chem(c_avg) - sig_ext (the lithiation path from
%   c_lo to c_hi, the delithiation path back), so they coincide.
%
%   Inputs (all optional)
%     params    parameter struct (default: default_parameters())
%     sig_ext   scalar or vector of scaled external stresses (default 0)
%     n_points  number of c_avg grid points (default 900)
%
%   Outputs (mu dimensionless, in units of R_g*T; Phi = -mu*R_g*T/F in mV;
%   rows correspond to sig_ext)
%     ocv.c_avg         1 x n           grid, c_avg in [c_lo, c_hi] = [0.005, 0.995]
%     ocv.mu_chem       1 x n           mu_chem(c_avg)
%     ocv.mu_single     n_sig x n       single-phase branch
%     ocv.single_state  1 x n           0 stable, 1 metastable, 2 unstable
%     ocv.mu_sep        n_sig x n       mixed-phase line (NaN outside [c_b1, c_b2])
%     ocv.mu_lith, ocv.mu_delith        n_sig x n   cycling paths on the grid
%     ocv.Phi_lith, ocv.Phi_delith      n_sig x n   the same as electrode potential
%     ocv.lith, ocv.delith              cycling paths in travel order, with fields
%         c (1 x m), s (1 x m, phase state), metastable (1 x m, logical),
%         mu (n_sig x m), Phi (n_sig x m)
%     ocv.phase_separation, ocv.sig_ext, ocv.c_b, ocv.c_sp, ocv.params
%
%   See also PLOT_OCV_HYSTERESIS, PLOT_GIBBS_MINIMISERS, PLOT_OCV_VOLTAGE.

if nargin < 1 || isempty(params),   params = default_parameters(); end
if nargin < 2 || isempty(sig_ext),  sig_ext = 0; end
if nargin < 3 || isempty(n_points), n_points = 900; end

params  = finalise_parameters(params);
params  = compute_critical_points(params);
sig_ext = sig_ext(:);
n_sig   = numel(sig_ext);
c_b     = params.c_b;
c_sp    = params.c_sp;
c_lo    = 0.005;
c_hi    = 0.995;
VT      = 1000*params.R_gT/params.F;        % thermal voltage [mV]
mu_of   = @(c, s, j) mu_hys(c, s, sig_ext(j), params)/params.R_gT;

%% Branches on a common grid
c_avg = linspace(c_lo, c_hi, n_points);
in_gap = c_avg >= c_b(1) & c_avg <= c_b(2);

single_state = zeros(1, n_points);
single_state((c_avg > c_b(1) & c_avg < c_sp(1)) | (c_avg > c_sp(2) & c_avg < c_b(2))) = 1;
single_state(c_avg >= c_sp(1) & c_avg <= c_sp(2)) = 2;

s_lith = zeros(1, n_points);
s_lith(c_avg > c_sp(1) & c_avg < c_b(2)) = 1;
s_delith = zeros(1, n_points);
s_delith(c_avg < c_sp(2) & c_avg > c_b(1)) = -1;

[mu_single, mu_sep, mu_lith, mu_delith] = deal(zeros(n_sig, n_points));
for j = 1:n_sig
    mu_single(j, :) = mu_of(c_avg, 0, j);
    mu_sep(j, :)    = mu_of(c_avg, 1, j);
    mu_lith(j, :)   = mu_of(c_avg, s_lith, j);
    mu_delith(j, :) = mu_of(c_avg, s_delith, j);
end
mu_sep(:, ~in_gap) = NaN;

%% Cycling paths in travel order, with the jumps at the spinodal points
m = max(20, round(n_points/4));             % points per segment
seg = @(a, b) linspace(a, b, m);

lith_c = [seg(c_lo, c_b(1)), seg(c_b(1), c_sp(1)), seg(c_sp(1), c_b(2)), seg(c_b(2), c_hi)];
lith_s = [zeros(1, 2*m), ones(1, m), zeros(1, m)];
lith_meta = [false(1, m), true(1, m), false(1, 2*m)];

delith_c = [seg(c_hi, c_b(2)), seg(c_b(2), c_sp(2)), seg(c_sp(2), c_b(1)), seg(c_b(1), c_lo)];
delith_s = [zeros(1, 2*m), -ones(1, m), zeros(1, m)];
delith_meta = [false(1, m), true(1, m), false(1, 2*m)];

if ~params.phase_separation                  % single phase throughout
    lith_c   = linspace(c_lo, c_hi, 4*m);
    delith_c = fliplr(lith_c);
    [lith_s, delith_s]       = deal(zeros(1, 4*m));
    [lith_meta, delith_meta] = deal(false(1, 4*m));
end

ocv.lith   = make_path(lith_c,   lith_s,   lith_meta,   n_sig, mu_of, VT);
ocv.delith = make_path(delith_c, delith_s, delith_meta, n_sig, mu_of, VT);

%% Pack
ocv.c_avg        = c_avg;
ocv.mu_chem      = mu_chem(c_avg);
ocv.mu_single    = mu_single;
ocv.single_state = single_state;
ocv.mu_sep       = mu_sep;
ocv.mu_lith      = mu_lith;
ocv.mu_delith    = mu_delith;
ocv.Phi_lith     = -VT*mu_lith;
ocv.Phi_delith   = -VT*mu_delith;
ocv.phase_separation = params.phase_separation;
ocv.sig_ext      = sig_ext;
ocv.c_b          = c_b;
ocv.c_sp         = c_sp;
ocv.params       = params;
end

% -------------------------------------------------------------------------
function path = make_path(c, s, metastable, n_sig, mu_of, VT)
path.c          = c;
path.s          = s;
path.metastable = metastable;
path.mu         = zeros(n_sig, numel(c));
for j = 1:n_sig
    path.mu(j, :) = mu_of(c, s, j);
end
path.Phi = -VT*path.mu;
end
