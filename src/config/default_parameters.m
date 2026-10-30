function params = default_parameters()
%DEFAULT_PARAMETERS Default parameter set for MechoCPM.
%
%   params = DEFAULT_PARAMETERS() returns a struct holding every model,
%   protocol and solver parameter. To run a different case, do not edit
%   this file -- override values instead:
%
%       sol    = run_simulation('Lambda', 1.32, 'sigma_y', 0.8, 'c_rate', 1/500);
%       params = set_parameters(default_parameters(), 'Lambda', 0.99);
%
%   Notation follows the paper (symbols in brackets). Stresses are scaled
%   by R_g*T/V_c, so the external stress enters the chemo-mechanical
%   potential with unit coefficient:  mu = mu_chem - sigma_kk/3.
%       Lambda  mechanical coupling parameter, Lambda = 2*G*V_c^2*c_m/(3*R_g*T)
%               (elastic misfit energy of a fully lithiated particle over R_g*T)
%       sigma_y scaled yield stress of the external shell, V_c*sigma_y/(R_g*T)
%       alpha   swelling magnitude, alpha = V_c*c_m/3 (displacement and strain
%               scale only; it does not enter the electrochemistry)
%       nu      Poisson ratio; nu_bar = nu/(1-2*nu) is derived (nu-bar)
%       c_b     binodal points [c_b1, c_b2]
%       c_sp    spinodal points [c_s1, c_s2]
%   Derived: A_int = Lambda*(1+3*nu_bar)/(3*(1+nu_bar)) (incompatibility
%   penalty) and Lambda_c = 3*(1+nu_bar)/(1+3*nu_bar), the coupling above
%   which phase separation is suppressed for the default mu_chem.
%
%   The previous scaling used gamma = 2*G*V_c/(R_g*T), with Lambda = gamma*alpha
%   and stresses scaled by 2*G. Old 'gamma' values raise an error with
%   the conversion; old sigma_y values must be multiplied by gamma.
%   Other older short names (alph, xstar, xspin, Far, n_total) are still
%   accepted by SET_PARAMETERS / RUN_SIMULATION.
%
%   See README.md for a full description of each parameter.

%% Physical constants
params.F   = 96500;             % Faraday constant [C/mol]
params.R_g = 8.314;             % universal gas constant [J/(mol K)]
params.T   = 300;               % temperature [K]

%% Host material
params.c_m = 2.7e5;            % maximum Li concentration [mol/m^3]
params.nu  = 0.25;             % Poisson ratio (nu_bar = nu/(1-2*nu) is derived)

%% Mechanics (main user-facing parameters)
params.Lambda  = 0;            % mechanical coupling Lambda (0 = no internal stress)
params.sigma_y = 0;            % scaled yield stress of the external shell, V_c*sigma_y/(R_g*T)
params.alpha   = 0.33;         % swelling magnitude alpha = V_c*c_m/3 (displacement scale only)
params.u0      = 1e-5;         % current scale setting the steepness of sig_ideal [A]
params.k       = 0.01;         % relaxation rate of sig_ext towards sig_ideal [1/s]

%% Reaction kinetics
params.I0 = 1;                 % exchange current density [A/m^2]

%% Particle population
params.n_groups          = 100;        % number of particle-size groups n
params.N_total           = 1e13;       % total number of particles N (scales Q and I only)
params.radius            = 70e-9;      % mean particle radius [m]
params.sigma_R           = 15e-9;      % half-width / spread of the size distribution [m]
params.size_distribution = 'uniform';  % 'uniform' or 'lognormal'

%% Cycling protocol
params.soc_waypoints = [0.01, 0.99, 0.01]; % SoC targets; one segment between each pair
params.c_rate        = 0.002;              % scalar, or one value per segment
params.rest_time     = 100;                % open-circuit rest before cycling [s]

%% Critical points: binodal c_b = [c_b1 c_b2], spinodal c_sp = [c_s1 c_s2]
params.compute_critical_points = true;     % false -> use c_b / c_sp below as given
params.c_b  = [0.0707, 0.9293];            % used only if compute_critical_points = false
params.c_sp = [0.2113, 0.7887];            % (values for mu_chem of the paper, Lambda = 0)

% Settings for FIND_BINODAL_POINTS
params.binodal_method        = 'analytical'; % 'analytical' (fsolve), 'gridsearch' or 'both'
params.binodal_grid_points   = 2000;         % grid points per axis ('gridsearch' / 'both')
params.binodal_tol           = 1e-3;         % max analytical-grid difference ('both')
params.binodal_initial_guess = [0.1, 0.9];   % fsolve starting point [c_c, c_s]
params.binodal_c_avg         = 0.51;         % c_avg used in the minimisation

%% Solver settings
params.RelTol               = 1e-6;
params.AbsTol               = 1e-8;
params.tau_I_default        = 2;     % current-ramp time constant [s]
params.tau_I_reversal       = 10;    % ramp time constant when the current changes sign [s]
params.max_events_per_block = 1e5;   % guard against event chattering
params.verbose              = true;

%% Derived quantities and validation
params.tau_I = params.tau_I_default;
params = finalise_parameters(params);
end
