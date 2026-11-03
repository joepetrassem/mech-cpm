function params = finalise_parameters(params)
%FINALISE_PARAMETERS Compute derived quantities and validate a parameter struct.
%
%   Derived fields:
%     nu_bar   = nu/(1 - 2*nu)
%     R_gT     = R_g*T
%     u0_used  = u0, or 0.05*min(c_rate)/3600 if u0 is empty (default) [1/s]
%     A_int    = Lambda*(1 + 3*nu_bar)/(3*(1 + nu_bar))   incompatibility penalty
%     Lambda_c = -min(mu_chem')*3*(1 + nu_bar)/(2*(1 + 3*nu_bar))
%                critical coupling: phase separation is possible only for
%                Lambda < Lambda_c, i.e. 2*A_int < -min(mu_chem'). Computed
%                from mu_chem (MIN_MU_CHEM_PRIME); 1.8 for the paper's
%                mu_chem (min mu_chem' = -2) and nu_bar = 0.5. Zero if
%                mu_chem is monotonic (no phase separation at all).
%     phase_separation = true if the critical points c_b and c_sp are finite.
%                It is false when mechanics suppresses phase separation, in
%                which case c_b = c_sp = [NaN NaN] and every particle stays
%                single phase (solid solution).

params.nu_bar   = params.nu/(1 - 2*params.nu);
params.R_gT     = params.R_g*params.T;
params.A_int    = params.Lambda*(1 + 3*params.nu_bar)/(3*(1 + params.nu_bar));
[~, dmu_min]    = min_mu_chem_prime();
params.Lambda_c = max(0, -dmu_min)*3*(1 + params.nu_bar)/(2*(1 + 3*params.nu_bar));

num = @(x, name, attrs) validateattributes(x, {'numeric'}, attrs, 'MechoCPM', name);

num(params.Lambda,   'Lambda',   {'scalar', 'real', 'finite', 'nonnegative'});
num(params.alpha,    'alpha',    {'scalar', 'real', 'finite', 'positive'});
num(params.sigma_y,  'sigma_y',  {'scalar', 'real', 'finite', 'nonnegative'});
num(params.nu,       'nu',       {'scalar', 'real', '>', -1, '<', 0.5});
num(params.k,        'k',        {'scalar', 'real', 'finite', 'nonnegative'});
num(params.I0,       'I0',       {'scalar', 'real', 'finite', 'positive'});
num(params.c_m,      'c_m',      {'scalar', 'real', 'finite', 'positive'});
num(params.F,        'F',        {'scalar', 'real', 'finite', 'positive'});
num(params.R_g,      'R_g',      {'scalar', 'real', 'finite', 'positive'});
num(params.T,        'T',        {'scalar', 'real', 'finite', 'positive'});
num(params.n_groups, 'n_groups', {'scalar', 'integer', 'positive'});
num(params.N_total,  'N_total',  {'scalar', 'real', 'finite', 'positive'});
num(params.radius,   'radius',   {'scalar', 'real', 'finite', 'positive'});
num(params.sigma_R,  'sigma_R',  {'scalar', 'real', 'finite', 'nonnegative'});
num(params.soc_waypoints, 'soc_waypoints', {'vector', 'real', '>', 0, '<', 1});
num(params.c_rate,   'c_rate',   {'vector', 'real', 'finite', 'positive'});
num(params.rest_time,'rest_time',{'scalar', 'real', 'finite', 'nonnegative'});
if isempty(params.u0)
    params.u0_used = 0.05*min(params.c_rate(:))/3600;    % 5% of the nominal swelling rate
else
    num(params.u0, 'u0', {'scalar', 'real', 'finite', 'positive'});
    params.u0_used = params.u0;
end
num(params.c_b,      'c_b',      {'numel', 2, 'real'});
num(params.c_sp,     'c_sp',     {'numel', 2, 'real'});
params.phase_separation = all(isfinite([params.c_b(:); params.c_sp(:)]));
if params.phase_separation
    num(params.c_b,  'c_b',  {'>', 0, '<', 1});
    num(params.c_sp, 'c_sp', {'>', 0, '<', 1});
elseif ~all(isnan([params.c_b(:); params.c_sp(:)]))
    error('MechoCPM:params', ...
        'c_b and c_sp must both be finite, or both [NaN NaN] for no phase separation.');
end
num(params.binodal_grid_points,   'binodal_grid_points',   {'scalar', 'integer', '>=', 10});
num(params.binodal_tol,           'binodal_tol',           {'scalar', 'real', 'positive'});
if ~isempty(params.binodal_initial_guess)
    num(params.binodal_initial_guess, 'binodal_initial_guess', {'numel', 2, 'real'});
end
if ~isempty(params.binodal_c_avg)
    num(params.binodal_c_avg, 'binodal_c_avg', {'scalar', 'real', '>', 0, '<', 1});
end

n_seg = numel(params.soc_waypoints) - 1;
if n_seg < 1
    error('MechoCPM:params', 'soc_waypoints needs at least two entries.');
end
if any(diff(params.soc_waypoints) == 0)
    error('MechoCPM:params', 'Consecutive soc_waypoints must differ.');
end
if ~(isscalar(params.c_rate) || numel(params.c_rate) == n_seg)
    error('MechoCPM:params', ...
        'c_rate must be a scalar or have one entry per segment (%d).', n_seg);
end
if ~any(strcmpi(params.size_distribution, {'uniform', 'lognormal'}))
    error('MechoCPM:params', 'size_distribution must be ''uniform'' or ''lognormal''.');
end
if ~isempty(params.max_step)
    num(params.max_step, 'max_step', {'scalar', 'real', 'positive'});
end
if ~any(strcmpi(params.solver, {'ode15s', 'ode15i'}))
    error('MechoCPM:params', 'solver must be ''ode15s'' or ''ode15i''.');
end
if ~any(strcmpi(params.binodal_method, {'analytical', 'gridsearch', 'both'}))
    error('MechoCPM:params', ...
        'binodal_method must be ''analytical'', ''gridsearch'' or ''both''.');
end
if strcmpi(params.size_distribution, 'uniform') && params.n_groups > 1 ...
        && params.radius - params.sigma_R <= 0
    error('MechoCPM:params', 'Uniform distribution needs radius > sigma_R.');
end
end
