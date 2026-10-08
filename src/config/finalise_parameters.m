function params = finalise_parameters(params)
%FINALISE_PARAMETERS Compute derived quantities and validate a parameter struct.
%
%   Derived fields: nu_bar = nu/(1 - 2*nu)  and  R_gT = R_g*T.

params.nu_bar = params.nu/(1 - 2*params.nu);
params.R_gT   = params.R_g*params.T;

num = @(x, name, attrs) validateattributes(x, {'numeric'}, attrs, 'MechoCPM', name);

num(params.alpha,    'alpha',    {'scalar', 'real', 'finite', 'positive'});
num(params.gamma,    'gamma',    {'scalar', 'real', 'finite', 'nonnegative'});
num(params.sigma_y,  'sigma_y',  {'scalar', 'real', 'finite', 'nonnegative'});
num(params.nu,       'nu',       {'scalar', 'real', '>', -1, '<', 0.5});
num(params.k,        'k',        {'scalar', 'real', 'finite', 'nonnegative'});
num(params.u0,       'u0',       {'scalar', 'real', 'finite', 'positive'});
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
num(params.c_b,      'c_b',      {'numel', 2, 'real', '>', 0, '<', 1});
num(params.c_sp,     'c_sp',     {'numel', 2, 'real', '>', 0, '<', 1});
num(params.binodal_grid_points,   'binodal_grid_points',   {'scalar', 'integer', '>=', 10});
num(params.binodal_tol,           'binodal_tol',           {'scalar', 'real', 'positive'});
num(params.binodal_initial_guess, 'binodal_initial_guess', {'numel', 2, 'real'});
num(params.binodal_c_avg,         'binodal_c_avg',         {'scalar', 'real', '>', 0, '<', 1});

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
if ~any(strcmpi(params.binodal_method, {'analytical', 'gridsearch', 'both'}))
    error('MechoCPM:params', ...
        'binodal_method must be ''analytical'', ''gridsearch'' or ''both''.');
end
if strcmpi(params.size_distribution, 'uniform') && params.n_groups > 1 ...
        && params.radius - params.sigma_R <= 0
    error('MechoCPM:params', 'Uniform distribution needs radius > sigma_R.');
end
end
