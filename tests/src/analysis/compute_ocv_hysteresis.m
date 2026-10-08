function ocv = compute_ocv_hysteresis(params, sig_ext, n_points)
%COMPUTE_OCV_HYSTERESIS Quasi-static single-particle potentials (Fig. 3a,b,d).
%
%   ocv = COMPUTE_OCV_HYSTERESIS(params, sig_ext, n_points) evaluates MU_HYS
%   on a grid of average concentrations c_avg for one or more external
%   stresses, without any time integration:
%     lithiation   : single phase until c_s1, then core-shell (s = +1, rich shell)
%                    until c_b2
%     delithiation : single phase until c_s2, then core-shell (s = -1, poor shell)
%                    until c_b1
%   The mixed-phase (global-minimum) curve of the paper is assembled from the
%   poor-shell branch for c_avg < 1/2 and the rich-shell branch for c_avg > 1/2.
%
%   Inputs (all optional)
%     params    parameter struct (default: default_parameters())
%     sig_ext   scalar or vector of external stresses (default 0)
%     n_points  number of c_avg points (default 900)
%
%   Outputs (dimensionless, units of R_g*T; rows correspond to sig_ext)
%     ocv.c_avg      1 x n_points         average concentration
%     ocv.mu_chem    1 x n_points         mu_chem(c_avg), the single-phase curve at gamma = 0
%     ocv.mu_lith    n_sig x n_points     lithiation branch
%     ocv.mu_delith  n_sig x n_points     delithiation branch
%     ocv.mu_mixed   n_sig x n_points     mixed-phase (equilibrium) curve
%     ocv.sig_ext, ocv.c_b, ocv.c_sp, ocv.params

if nargin < 1 || isempty(params),   params = default_parameters(); end
if nargin < 2 || isempty(sig_ext),  sig_ext = 0; end
if nargin < 3 || isempty(n_points), n_points = 900; end

params = finalise_parameters(params);
params = compute_critical_points(params);

sig_ext = sig_ext(:);
c_avg   = linspace(0.005, 0.995, n_points);

s_lith = zeros(size(c_avg));
s_lith(c_avg > params.c_sp(1) & c_avg < params.c_b(2)) = 1;
s_delith = zeros(size(c_avg));
s_delith(c_avg < params.c_sp(2) & c_avg > params.c_b(1)) = -1;

n_sig = numel(sig_ext);
mu_lith   = zeros(n_sig, n_points);
mu_delith = zeros(n_sig, n_points);
for j = 1:n_sig
    mu_lith(j, :)   = mu_hys(c_avg, s_lith,   sig_ext(j), params)/params.R_gT;
    mu_delith(j, :) = mu_hys(c_avg, s_delith, sig_ext(j), params)/params.R_gT;
end

half = floor(n_points/2);
mu_mixed = [mu_delith(:, 1:half), mu_lith(:, half+1:end)];

ocv.c_avg     = c_avg;
ocv.mu_chem   = mu_chem(c_avg);
ocv.mu_lith   = mu_lith;
ocv.mu_delith = mu_delith;
ocv.mu_mixed  = mu_mixed;
ocv.sig_ext   = sig_ext;
ocv.c_b       = params.c_b;
ocv.c_sp      = params.c_sp;
ocv.params    = params;
end
