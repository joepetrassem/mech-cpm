function [c_b, info] = find_binodal_points(params, varargin)
%find_binodal_points Binodal (coexisting core and shell) concentrations.
%
%   c_b = FIND_BINODAL_POINTS(params) returns the binodal concentrations
%   c_b = [c_b1, c_b2] (dimensionless, ascending) of a phase-separated
%   particle: the core and shell concentrations that minimise its Gibbs free
%   energy. By default they are found by solving the stationarity conditions
%   of the Gibbs free energy with FSOLVE (the analytical solution).
%
%   c_b = find_binodal_points(params, Name, Value) sets options:
%     'Method'        'analytical' (default) - solve the stationarity
%                     conditions with FSOLVE.
%                     'gridsearch' - minimise the Gibbs free energy by brute
%                     force over a grid of core and shell concentrations.
%                     'both' - compute both and compare them; the grid-search
%                     result is returned (with a warning) if they differ by
%                     more than 'Tolerance'.
%     'GridPoints'    Grid points per concentration axis (default 2000).
%     'Tolerance'     Maximum allowed difference between the two methods
%                     (default 1e-3). Only used with 'both'.
%     'InitialGuess'  Starting point for FSOLVE (default [0.1 0.9]).
%     'c_avg'         Average particle concentration used to evaluate the
%                     core volume fraction (default 0.51). In the
%                     sharp-interface, linear-elastic model the binodal
%                     points do not depend on it, provided it lies between
%                     them.
%     'Verbose'       Print a summary and timing (default false).
%
%   [c_b, info] = find_binodal_points(...) also returns a structure with fields
%     method       Method used
%     analytical   FSOLVE solution ([] if not computed)
%     exitflag     FSOLVE exit flag ([] if not computed)
%     grid         Grid-search solution ([] if not computed)
%     difference   analytical - grid ([] unless Method is 'both')
%
%   The external stress is set to zero: the binodal points are independent
%   of it in this model.
%
%   Requires the Optimization Toolbox (FSOLVE) for 'analytical' and 'both'.
%
%   The analytical method solves Eqs. (59)-(62) of the paper for the core and
%   shell concentrations (c_c, c_s); the sorted pair is [c_b1, c_b2].
%
%   See also FIND_SPINODAL_POINTS, COMPUTE_CRITICAL_POINTS, F_CHEM, MU_CHEM.

% ---------------------------------------------------------------- options --
p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'params', @isstruct);
addParameter(p, 'Method', 'analytical', @(x) ischar(x) || isstring(x));
addParameter(p, 'GridPoints', 2000, @(x) isnumeric(x) && isscalar(x) && x >= 10 && x == round(x));
addParameter(p, 'Tolerance', 1e-3, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'InitialGuess', [0.1 0.9], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'c_avg', 0.51, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'Verbose', false, @(x) islogical(x) || isnumeric(x));
parse(p, params, varargin{:});
opt = p.Results;
method = validatestring(opt.Method, {'analytical', 'gridsearch', 'both'}, mfilename, 'Method');

% ------------------------------------------------- model coefficients --
mech.gamma = params.gamma;          % mechanical coupling parameter gamma
mech.alpha = params.alpha;          % swelling parameter alpha
mech.cons  = 1 + 3*params.nu_bar;   % 1 + 3*nu_bar
mech.c_avg = opt.c_avg;

info = struct('method', method, 'analytical', [], 'exitflag', [], ...
              'grid', [], 'difference', []);
timer = tic;

% ------------------------------------------------------------- solve --
switch method
    case 'analytical'
        [c_b, info.exitflag] = solve_analytical(mech, opt.InitialGuess);
        info.analytical = c_b;

    case 'gridsearch'
        c_b = solve_grid(mech, opt.GridPoints);
        info.grid = c_b;

    case 'both'
        [c_ana, info.exitflag] = solve_analytical(mech, opt.InitialGuess);
        c_grid = solve_grid(mech, opt.GridPoints);
        info.analytical = c_ana;
        info.grid = c_grid;
        info.difference = c_ana - c_grid;
        if any(abs(info.difference) > opt.Tolerance)
            warning('find_binodal_points:methodsDisagree', ...
                ['Analytical [%.4f %.4f] and grid-search [%.4f %.4f] binodal points ', ...
                 'differ by more than %.1e; returning the grid-search result.'], ...
                c_ana, c_grid, opt.Tolerance);
            c_b = c_grid;
        else
            c_b = c_ana;
        end
end

if opt.Verbose
    fprintf('Binodal points (%s): c_b1 = %.6f, c_b2 = %.6f  [%.2f s]\n', ...
        method, c_b(1), c_b(2), toc(timer));
    if ~isempty(info.difference)
        fprintf('  analytical - grid = [%.2e %.2e]\n', info.difference);
    end
end
end

% =========================================================================
function [c_b, exitflag] = solve_analytical(mech, guess)
% Solve the two stationarity conditions of the Gibbs free energy for the
% core concentration c(1) and shell concentration c(2).
g = mech.gamma; a = mech.alpha; k = mech.cons; c_avg = mech.c_avg;

psi    = @(c) (c_avg - c(2)) ./ (c(1) - c(2));                   % core volume fraction
lambda = @(c) psi(c)*mu_chem(c(1)) + (1 - psi(c))*mu_chem(c(2)); % Lagrange multiplier

residual = @(c) [ ...
    mu_chem(c(1)) - 2*g*a*k*(c(2) - c(1))*(1 - psi(c))/(k + 2) - lambda(c); ...
    f_chem(c(1)) - f_chem(c(2)) + g*a*k*(c(2) - c(1))^2*(1 - 2*psi(c))/(k + 2) ...
        + lambda(c)*(c(2) - c(1))];

options = optimoptions('fsolve', 'Display', 'off', ...
    'Algorithm', 'trust-region-dogleg', 'OptimalityTolerance', 1e-11);
[c, ~, exitflag] = fsolve(residual, guess(:).', options);
c_b = sort(c, 'ascend');

if exitflag <= 0 || any(c_b <= 0 | c_b >= 1) || diff(c_b) < 1e-5
    error('find_binodal_points:analyticalFailed', ...
        ['Analytical solution failed (exitflag %d, c = [%.4g %.4g]). ', ...
         'Try another ''InitialGuess'' or ''Method'', ''gridsearch''.'], ...
        exitflag, c_b(1), c_b(2));
end
end

% =========================================================================
function c_b = solve_grid(mech, n)
% Minimise the Gibbs free energy over an n-by-n grid of core and shell
% concentrations. The grid excludes 0 and 1, where f_chem is singular.
% Loops over shell concentrations and vectorises over core concentrations,
% which keeps memory use at O(n).
g = mech.gamma; a = mech.alpha; k = mech.cons; c_avg = mech.c_avg;
sig = 0;                                   % external stress (no effect on c_b)

c  = linspace(0, 1, n + 2);  c = c(2:end-1);
fc = arrayfun(@f_chem, c);                 % chemical free energy on the grid

G_best = Inf; c_b = [NaN NaN];
for i = 1:n
    c_core  = c;                           % vector of core concentrations
    c_shell = c(i);                        % one shell concentration
    D   = c_shell - c_core;
    psi = (c_avg - c_shell) ./ (c_core - c_shell);

    w_core  = 3*(sig + a*k*D.*2.*(1 - psi)/(k + 2)) .* (sig/(a*k) + D.*2.*(1 - psi)/(k + 2));
    w_shell = 3*(sig^2/(a*k) + 2*sig*D.*2.*(-psi)/(k + 2) ...
                 + 2*(D/(k + 2)).^2*a*k.*psi.*(2*psi + k));

    g_core  = fc       + g*0.5*w_core/3;
    g_shell = fc(i)    + g*0.5*w_shell/3;

    % The Lagrange term lambda*(psi*c_core + (1-psi)*c_shell - c_avg) is
    % identically zero by the definition of psi, so it is omitted.
    G = psi.*g_core + (1 - psi).*g_shell;
    G(~isfinite(psi) | psi < 0 | psi > 1) = NaN;

    [G_row, j] = min(G);                   % min ignores NaN
    if G_row < G_best
        G_best = G_row;
        c_b = [c_core(j), c_shell];
    end
end

if any(isnan(c_b))
    error('find_binodal_points:gridFailed', ...
        'Grid search found no admissible phase-separated state.');
end
c_b = sort(c_b, 'ascend');
end
