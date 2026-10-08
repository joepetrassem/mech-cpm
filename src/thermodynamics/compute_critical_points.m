function params = compute_critical_points(params)
%COMPUTE_CRITICAL_POINTS Fill params.c_b (binodal) and params.c_sp (spinodal).
%
%   If params.compute_critical_points is true (default), the binodal points
%   c_b = [c_b1, c_b2] come from FIND_BINODAL_POINTS (method chosen by
%   params.binodal_method: 'analytical', 'gridsearch' or 'both') and the
%   spinodal points c_sp = [c_s1, c_s2] from FIND_SPINODAL_POINTS. Otherwise
%   the values already in params are used unchanged.

if params.compute_critical_points
    c_sp = find_spinodal_points(params);
    if any(isnan(c_sp))
        error('MechoCPM:noSpinodal', ...
            ['No spinodal region for gamma = %g, alpha = %g: mechanics suppresses ' ...
             'phase separation, so there are no binodal points either. Reduce gamma.'], ...
            params.gamma, params.alpha);
    end
    c_b = find_binodal_points(params, ...
        'Method',       params.binodal_method, ...
        'GridPoints',   params.binodal_grid_points, ...
        'Tolerance',    params.binodal_tol, ...
        'InitialGuess', params.binodal_initial_guess, ...
        'c_avg',        params.binodal_c_avg);
    if diff(sort(c_b)) < 1e-5
        error('MechoCPM:binodalPoints', ...
            'Binodal points are (nearly) coincident: [%g %g].', c_b);
    end
    params.c_b  = c_b;
    params.c_sp = c_sp;
end
params.c_b  = sort(params.c_b(:)');
params.c_sp = sort(params.c_sp(:)');

if ~(params.c_b(1) < params.c_sp(1) && params.c_sp(2) < params.c_b(2))
    warning('MechoCPM:criticalOrder', ...
        'Expected c_b1 < c_s1 < c_s2 < c_b2; got c_b = [%.4f %.4f], c_sp = [%.4f %.4f].', ...
        params.c_b, params.c_sp);
end
end
