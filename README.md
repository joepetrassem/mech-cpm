# MechoCPM

**Chemo-mechanical model of phase-changing electrode particles**

MATLAB code accompanying *Coupled chemo-mechanical framework for phase-changing electrodes* (J. Petrassem de Sousa, W. Clarke, A. Galvis, S. Sahu and J. M. Foster, Journal of The Electrochemical Society, 2026).

The code simulates galvanostatic cycling of an ensemble of electrode particles with a distribution of sizes. Each particle is either single-phase or phase-separated into a core and a shell. Transitions are triggered when the particle's average concentration `c_avg` crosses a spinodal point (single phase -> core-shell) or a binodal point (core-shell -> single phase). Mechanics enters in two independent ways, following the scaling of the paper. The mechanical coupling `Lambda`, the ratio of the elastic misfit energy of a fully lithiated particle to the thermal energy, sets the internal incompatibility stress between core and shell: it shifts the spinodal and binodal points and tilts the mixed-phase branch of the potential. An external radial stress `sigma_ext`, which relaxes towards a plastic yield stress `sigma_y`, translates the potential rigidly. Stresses are measured in units of `R_g*T/V_c` (about 0.3 GPa), so the external stress enters the potential with unit coefficient.

Symbols in the code follow the paper (see [Notation](#notation)); the source comments name the paper sections they implement.

## Contents

- [Quick start](#quick-start)
- [Testing that everything works](#testing-that-everything-works)
- [Choosing parameters](#choosing-parameters)
- [Binodal points: analytical vs grid search](#binodal-points-analytical-vs-grid-search)
- [Notation](#notation)
- [Model summary](#model-summary)
- [External medium](#external-medium)
- [Time integration](#time-integration)
- [Using your own chemical potential](#using-your-own-chemical-potential)
- [Repository layout](#repository-layout)
- [Examples](#examples)
- [Output](#output)
- [Publishing this repository on GitHub](#publishing-this-repository-on-github)

## Quick start

Requirements: MATLAB R2019b or newer. The Optimization Toolbox (`fsolve`) is needed for the default analytical binodal solver; `ode15s`, `ode15i` and `fzero` are in base MATLAB, and the grid-search binodal method needs no toolbox.

```matlab
cd path/to/MechoCPM          % the folder containing this README
setup_paths();               % adds src/ to the MATLAB path

sol = run_simulation('Lambda', 1.32, 'sigma_y', 0.8, 'n_groups', 20, 'c_rate', 0.01);

plot_voltage_soc(sol);          % electrode potential [mV] against SoC
plot_voltage_time(sol);         % ... and against time
voltage_summary(sol);           % prints Phi on each branch and the hysteresis gap
plot_phase_fraction(sol);
```

## Testing that everything works

Run these from the repository root in MATLAB (the folder that contains `setup_paths.m`).

**1. Unit and regression tests (about 20 s).**

```matlab
setup_paths();
results = runtests('tests');
disp(table(results))
```

All 54 tests should pass (`test_examples` runs examples 1-6 and 8 without figures on screen, so allow a few minutes). They check, among other things, that

- `mu_chem`, `mu_chem_prime` and `f_chem` are mutually consistent and equal the paper's `mu_chem`;
- the spinodal points satisfy the coherent-spinodal condition `mu_chem'(c) + 2*A_int = 0`, exist only for `Lambda < Lambda_c` (1.8 for the paper's `mu_chem`), match the closed form of the regular solution, and start-up inside the spinodal range is rejected; and the binodal points satisfy the common-tangent construction at `Lambda = 0`;
- the analytical and grid-search binodal methods agree, and the binodal points do not depend on `c_avg`;
- `mu_hys` equals the paper's `mu_single` and `mu_sep`, the core and shell potentials are equal at equilibrium, and the external stress shifts the potential by exactly `-sig_ext` for any `Lambda`;
- old `gamma` inputs are rejected with the conversion to `Lambda`, and examples 1-6 run and produce figures and logs that never mention `gamma`;
- the plastic-shell external stress starts at zero, opposes the surface motion (`-sigma_y` by the end of lithiation, `+sigma_y` by the end of delithiation), stays continuous across blocks, unloads elastically over `Delta c = 2*sigma_y/(Lambda*(1+3*nu_bar))` after a reversal, and is zero for `Lambda = 0`; the implicit stress law is solved exactly, and the DAE Jacobian supplied to `ode15i` matches finite differences;
- the elastic-shell external stress follows `sig_ext = -K_eff*(c_avg - c_form)` along a whole cycle (both solvers), shifts the potential on both branches alike without adding a gap, and gives the expected `K_eff` limits;
- the reduced ODE solved with `ode15s` (the default) satisfies every equation of the DAE, and gives the same electrode potential, events and stresses as `ode15i`;
- once the stress has switched, it widens the simulated voltage gap by exactly `2*sigma_y*R_g*T/F`, and the quasi-static equilibrium potential is `-mu*R_g*T/F`;
- the single-particle cycling paths (Fig. 3c,d) jump onto the mixed-phase line at the spinodal points and rejoin the single-phase branch continuously at the binodal points, single-phase states are classified correctly as stable, metastable or unstable, and the rich-shell and poor-shell mixed-phase branches coincide; the arrowheads on those plots lie on the paths, tangent to them on screen and pointing in the direction of travel;
- above `Lambda_c` the particle is a solid solution: the quasi-static paths coincide on `mu_chem - sig_ext` for any `Lambda`, and a simulation runs with no phase-change events while still conserving charge;
- the DAE initial conditions are consistent and a short cycle conserves charge;
- **`test_regression_against_original_code`**: the electrode potential, binodal and spinodal points reproduce results saved from the original, pre-reorganisation code (`tests/data/original_code_reference.mat`) for the two reference parameter sets without external stress, to within solver tolerance (RMS difference below 10 microvolts). The two references with external stress are not used, because the original code drove the stress by the global current and reset it at every block. The reference was made with the previous scaling and is converted with `Lambda = gamma*alpha`, `sigma_y = gamma*sigma_y_old`, which checks that the change of scaling leaves the model unchanged.

**2. Examples.** Run them one at a time (each opens figures) or all at once:

```matlab
cd examples
run_all_examples();        % examples 1-6, about 40 seconds
run_all_examples(true);    % also example 7 (paper Figs. 4-7), another 1-2 minutes
```

What to expect from each example:

| Example | Runtime | What you should see |
|---|---|---|
| `example_01_single_cycle` | seconds | 6 figures: electrode potential against SoC (a plateau with small steps) and against time, particle potentials, external stress switching between about -0.6 and +0.6, and the fraction of phase-separated particles rising and falling. Prints the potential on each branch at SoC = 0.25, 0.5, 0.75 and the gap |
| `example_02_Lambda_sweep` | seconds | Electrode potential against SoC and time, and phase fractions, for `Lambda` = 0, 0.33, 0.66, 1.32; the plateau shrinks and tilts as `Lambda` grows. Prints the potential and gap for each run |
| `example_03_sigma_y_sweep` | seconds | Electrode potential against SoC and time for `sigma_y` = 0, 0.1, 0.2, 0.4 with `Lambda = 1.32`: larger `sigma_y` widens the gap between lithiation and delithiation. Prints each switching window `2*sigma_y/(Lambda*(1+3*nu_bar))` (at most 0.24) and the simulated increase of the gap at SoC = 0.5 next to the prediction `2*sigma_y*R_g*T/F` (about 52 mV per unit); the two agree because every particle has switched by then |
| `example_04_crate_and_particles` | seconds | Electrode potential against SoC and time, with printed summaries, for a C-rate sweep, a number-of-groups sweep and per-segment C-rates |
| `example_05_ocv_hysteresis` | 1 s | Paper Fig. 3a,b: single-phase (dark blue) and mixed-phase (blue) solution branches for `Lambda = 0` (flat mixed-phase line) and `Lambda = 0.99` (sloped line). Fig. 3c,d: lithiation and delithiation paths for `Lambda = 0` and `Lambda = 0.99`, each with `sig_ext = -0.6` (red), 0 (black), 0.6 (blue): single phase up to the spinodal point (metastable beyond the binodal point), a vertical jump onto the mixed-phase line, phase separation to the far binodal point, the unstable single-phase states dotted, and arrowheads tangent to the paths showing the direction of charge. Also `Lambda = 2 > Lambda_c`, where there is no phase separation and both paths follow `mu_chem - sig_ext`. Prints `c_b = [0.0707 0.9293]`, `c_sp = [0.2113 0.7887]` and jumps in mu of -0.415 / +0.415 for `Lambda = 0`; `c_b = [0.1570 0.8430]`, `c_sp = [0.2857 0.7143]` and jumps of -0.134 / +0.134 for `Lambda = 0.99`; and the equilibrium potentials at c_avg = 0.5, which shift by about +-15.5 mV for `sig_ext = -+0.6` |
| `example_06_critical_points_vs_Lambda` | seconds | Paper Fig. 2: spinodal and binodal points against `Lambda`, with `Lambda_c = 1.8` marked; prints the largest analytical vs grid-search difference (below 1e-3) |
| `example_07_external_stress_comparison` | several minutes | Paper Figs. 4/6 (filling ratio of each particle group), Fig. 7 (electrode potential against SoC with and without external stress, `Lambda = 1.32`, `sigma_y = 0.8`), the electrode potential against time and the external stress of each particle. With these values the stress takes `Delta c` of about 0.48 to switch after the reversal, so the full gap increase of about 41 mV (`2*0.8*R_g*T/F`) is only reached on the second half of delithiation; the script prints the gap increase at SoC = 0.25, 0.5 and 0.75 |
| `example_08_elastic_shell` | about a minute | Electrode potential against SoC and time for an elastic shell with `k_e` = 0, 2, 10 and `Inf` (rigid), `Lambda = 0.99`: the curves tilt more steeply as `k_e` grows, with no extra gap between lithiation and delithiation. Prints `K_eff` for each run and confirms `sig_ext = -K_eff*(c_avg - c_form)` for every particle; plots the stress against `c_avg` and the quasi-static curves with and without the shell |

**3. A quick check by hand.** With no mechanics the binodal points are the classical regular-solution values:

```matlab
setup_paths();
p = compute_critical_points(default_parameters());
p.c_b      % 0.0707  0.9293
p.c_sp     % 0.2113  0.7887
```

## Choosing parameters

All parameters live in a single struct created by `default_parameters()`. Override any of them by name, either directly in `run_simulation` or with `set_parameters`:

```matlab
sol = run_simulation('Lambda', 0.66, 'c_rate', 0.05);          % name-value pairs

params = default_parameters();
params = set_parameters(params, 'Lambda', 0.99, 'sigma_y', 0.6); % build a struct...
sol    = run_simulation(params, 'n_groups', 50);               % ...and override further
```

Names are case-insensitive, misspelt names raise an error listing the valid ones, and values are validated before anything is solved. Names from earlier versions of the code (`alph`, `sig_y`, `Far`, `xstar`, `xspin`, `n_total`) are still accepted. `gamma` (and `gamm`) from the previous scaling raise an error that explains the conversion below.

### Main parameters

| Name | Default | Paper symbol | Description |
|---|---|---|---|
| `Lambda` | 0 | Lambda | Mechanical coupling, `2*G*V_c^2*c_m/(3*R_g*T)`: elastic misfit energy of a fully lithiated particle over `R_g*T` (0 = no internal stress; phase separation needs `Lambda < Lambda_c = 1.8`) |
| `external_medium` | `'plastic'` | | External medium around each particle: `'plastic'` (rapidly yielding plastic shell), `'elastic'` (elastic shell) or `'none'`. See [External medium](#external-medium) |
| `sigma_y` | 0 | sigma_y^* | Plastic shell: scaled yield stress, `V_c*sigma_y/(R_g*T)` |
| `k_e` | 0 | k_e | Elastic shell: scaled stiffness, the same for all particles (`Inf` = rigid medium) |
| `alpha` | 0.33 | alpha | Swelling magnitude `V_c*c_m/3`; displacement and strain scale only, does not affect the electrochemistry |
| `nu` | 0.25 | nu | Poisson ratio; `nu_bar = nu/(1-2*nu)` (= 0.5) is derived |
| `n_groups` | 100 | n | Number of particle-size groups |
| `c_rate` | 0.002 | | C-rate; a scalar, or one value per protocol segment |
| `soc_waypoints` | `[0.01 0.99 0.01]` | | SoC targets; the protocol has one segment between each pair |

### Other parameters

| Name | Default | Description |
|---|---|---|
| `radius`, `sigma_R` | 70 nm, 15 nm | Mean radius R-bar and spread of the size distribution |
| `size_distribution` | `'uniform'` | `'uniform'` or `'lognormal'` |
| `N_total` | 1e13 | Total number of particles N (only rescales capacity and current) |
| `rest_time` | 100 s | Open-circuit rest before cycling |
| `c_form` | `[]` (auto) | Elastic shell: filling fraction at which the shell was formed stress-free; empty means the initial filling `soc_waypoints(1)` |
| `k` | 0.01 1/s | Plastic shell: relaxation rate of the external stress |
| `u0` | `[]` (auto) | Plastic shell: reference surface velocity in `sigma_ideal` [1/s]. Empty: `0.05*min(c_rate)/3600`, 5% of the nominal swelling rate (stored as `u0_used`). The shell acts as rigid-plastic when `u0` is much smaller than the swelling rate and `Lambda*(1+3*nu_bar)*u0` much smaller than `k*sigma_y`; `run_simulation` warns if not. Much smaller values make the switch very sharp, which slows the integration (and can stop `ode15i`) |
| `I0` | 1 A/m^2 | Exchange current density |
| `c_m` | 2.7e5 mol/m^3 | Maximum concentration |
| `F`, `R_g`, `T` | 96500, 8.314, 300 | Faraday constant, gas constant, temperature |
| `compute_critical_points` | `true` | If `false`, `c_b` and `c_sp` are used as given |
| `binodal_method` | `'analytical'` | `'analytical'`, `'gridsearch'` or `'both'` (see below) |
| `binodal_grid_points`, `binodal_tol` | 2000, 1e-3 | Grid resolution and agreement tolerance for the grid search |
| `binodal_initial_guess`, `binodal_c_avg` | `[]`, `[]` (auto) | `fsolve` starting point and `c_avg` used in the minimisation; empty means a guess built from the spinodal points and `c_avg = mean(c_sp)` |
| `solver` | `'ode15s'` | `'ode15s'`: reduced ODE in `(c_avg, sig_ext)` (default); `'ode15i'`: fully implicit DAE. See [Time integration](#time-integration) |
| `RelTol`, `AbsTol` | 1e-6, 1e-8 | Integration tolerances |
| `max_step` | `[]` (auto) | Largest `ode15s` step [s]; empty means 1/1000 of each block, so the output resolves the voltage curve |
| `tau_I_default`, `tau_I_reversal` | 2 s, 10 s | Time constants of the current ramp between blocks |
| `verbose` | `true` | Print progress |

Derived (read-only) fields: `nu_bar = nu/(1-2*nu)`, `R_gT = R_g*T`, the incompatibility penalty `A_int = Lambda*(1+3*nu_bar)/(3*(1+nu_bar))` and the critical coupling `Lambda_c = 3*(1+nu_bar)/(1+3*nu_bar)`.

### Converting from the previous scaling

Earlier versions used `gamma = 2*G*V_c/(R_g*T)` and measured stresses in units of `2*G`. The model is unchanged; only the parameters are regrouped:

| Previous | Now | Paper values |
|---|---|---|
| `gamma` | `Lambda = gamma*alpha` | `gamma = 3, 4` -> `Lambda = 0.99, 1.32` |
| `sigma_y` | `sigma_y = gamma*sigma_y_old` | `gamma = 4, sigma_y = 0.2` -> `sigma_y = 0.8` |
| `mu = mu_chem - gamma*sig_ext` | `mu = mu_chem - sig_ext` | |

With the new scaling the external stress acts even when `Lambda = 0`, which the previous form (where `gamma = 0` also removed the external stress) could not express.

### Parameter sweeps

```matlab
[sols, labels] = run_parameter_sweep('Lambda', [0 0.33 0.66 1.32], 'n_groups', 10);
plot_voltage_soc(sols, labels);
voltage_summary(sols, labels);
```

## Binodal points: analytical vs grid search

The binodal points `c_b = [c_b1, c_b2]` are the core and shell concentrations that minimise the Gibbs free energy of a phase-separated particle. `find_binodal_points` offers two independent ways to compute them:

| Method | What it does | Needs |
|---|---|---|
| `'analytical'` (default) | Solves the reduced stationarity conditions (paper, Appendix D) with `fsolve` | Optimization Toolbox |
| `'gridsearch'` | Minimises `G = psi*g_chem(c_c) + (1-psi)*g_chem(c_s) + A_int*psi*(1-psi)*(c_s-c_c)^2` by brute force over a grid of core and shell concentrations | nothing extra |
| `'both'` | Computes both and compares them; warns and returns the grid result if they differ by more than the tolerance | Optimization Toolbox |

```matlab
p = set_parameters(default_parameters(), 'Lambda', 0.99);
c_b = find_binodal_points(p);                                   % analytical
c_b = find_binodal_points(p, 'Method', 'gridsearch');
[c_b, info] = find_binodal_points(p, 'Method', 'both', 'Verbose', true);

sol = run_simulation('Lambda', 0.99, 'binodal_method', 'both'); % choose it for a simulation
```

The spinodal points `c_sp = [c_s1, c_s2]` come from the coherent-spinodal condition `mu_chem'(c) + 2*A_int = 0` in `find_spinodal_points`. For `Lambda >= Lambda_c` (computed from `mu_chem`; 1.8 for the paper's `mu_chem` and `nu = 0.25`) no spinodal points exist: phase separation is suppressed. The critical points are then returned as `NaN`, `params.phase_separation` is `false`, and every particle lithiates as a solid solution with `mu = mu_chem(c_avg) - sig_ext`. Lithiation and delithiation coincide (apart from any hysteretic external stress), the result no longer depends on `Lambda`, and the decreasing part of `mu_chem` is stable. `run_simulation` runs without phase-change events. The external stress does not affect either set of critical points.

## Notation

| Paper | Code | Meaning |
|---|---|---|
| Lambda, alpha | `Lambda`, `alpha` | Mechanical coupling parameter, swelling magnitude |
| A_int | `A_int` | Incompatibility penalty, `Lambda*(1+3nu_bar)/(3(1+nu_bar))` |
| nu, nu-bar | `nu`, `nu_bar` | Poisson ratio, `nu/(1-2nu)` |
| sigma_y^* | `sigma_y` | Scaled yield stress of the external shell |
| sigma_rr^ext | `sig_ext` | External radial stress at the particle surface, in units of `R_g*T/V_c` |
| c_avg | `c_av` (dimensional, mol/m^3), `theta = c_av/c_m` or `c_avg` (dimensionless) | Average concentration |
| c_c, c_s | `c_c`, `c_s` | Core and shell concentrations |
| psi | `psi` | Core volume fraction |
| c_b1, c_b2 | `c_b(1)`, `c_b(2)` | Binodal points |
| c_s1, c_s2 | `c_sp(1)`, `c_sp(2)` | Spinodal points |
| mu_chem(c), g_chem(c) | `mu_chem`, `f_chem` | Homogeneous chemical potential and free energy |
| mu | `mu_hys` | Surface chemo-mechanical potential (`mu_single`, `mu_sep`) |
| Phi | `Phi` | Electrode potential |
| j_tr | `jtr` | Interfacial transfer current density |
| N^(i), A^(i), V^(i) | `N_part`, `Apart`, `Vol` | Particles, surface area and volume per group |
| k, u0, I0, c_m | `k`, `u0`, `I0`, `c_m` | Stress relaxation rate, reference surface velocity, exchange current density, maximum concentration |
| u_dot | `surface_velocity` | Surface velocity of a particle |

The phase state `s` of a particle is `0` (single phase), `+1` (lithium-rich shell, core `c_b1`, shell `c_b2`) or `-1` (lithium-poor shell, core `c_b2`, shell `c_b1`).

## Model summary

For size group `i` (radius `R_i`) with average concentration `c_avg^(i)`, phase state `s_i`, external stress `sig_ext^(i)` and a common electrode potential `Phi`:

```
dc_avg^(i)/dt  = -j_tr^(i) A^(i) / (F V^(i))
eta^(i)        = (F Phi + mu^(i)) / (R_g T)
j_tr^(i)       = 2 I0 sqrt(c_surf (c_m - c_surf)) / c_m * sinh(eta/2)
I(t)           = sum_i N^(i) A^(i) j_tr^(i)
external stress: plastic shell, elastic shell or none (see External medium)
```

`mu^(i)` is the surface chemo-mechanical potential from the single-particle model (`mu_hys`): in units of `R_g*T`,

```
mu_single = mu_chem(c_avg) - sig_ext
mu_sep    = mu_chem(c_s) - 2 A_int (c_avg - c_s) - sig_ext,   A_int = Lambda (1+3 nu_bar) / (3 (1+nu_bar))
```

for a single-phase and a core-shell particle (shell concentration `c_s`). `Lambda` tilts the mixed-phase branch; the external stress translates both branches by `-sig_ext`, so a hysteretic stress of `+-sigma_y` widens the voltage gap by `2*sigma_y*R_g*T/F`. The current is ramped smoothly between blocks (`ramp_current`).

### External medium

The external radial stress `sig_ext` enters the single-particle problem only as a boundary condition, so different surrounding media can be used without changing the minimisation. Choose one with `external_medium`:

**`'plastic'` (default): rapidly yielding plastic shell.** Parameters `sigma_y`, `k`, `u0`.

```
d sig_ext^(i)/dt = -k (sig_ext^(i) - sig_ideal^(i))
sig_ideal^(i)    = -sigma_y tanh(u_dot^(i)/u0)
u_dot^(i)        = (d sig_ext^(i)/dt) / (Lambda (1+3 nu_bar)) + (dc_avg^(i)/dt) / c_m
```

The shell (`sig_ideal`, `surface_velocity`). It responds to the velocity `u_dot` of that particle's own surface: while `|sig_ext| < sigma_y` the shell holds the surface still and the stress builds at a rate set by the particle stiffness; once it reaches the yield value the shell flows and the stress stays at `-sigma_y` (expanding) or `+sigma_y` (contracting). The law is implicit in `d sig_ext/dt`; `solve_stress_rate` solves it for the stress rate (see [Time integration](#time-integration)). In an ensemble, a particle that pauses or gives back lithium while others react unloads partly and then reloads, so during lithiation the stresses are not all exactly `-sigma_y`. Each particle starts unstressed and the stress is continuous across current blocks, so after a reversal it unloads elastically over `Delta c = 2*sigma_y/(Lambda*(1+3*nu_bar))` of filling fraction. For `Lambda = 0` the particle never loads the shell and `sig_ext = 0`. To model a different surrounding medium, change `sig_ideal` (both its value and its derivative) and, if needed, the relaxation law in `solve_stress_rate`, `dae_residual` and `dae_jacobian`.

**`'elastic'`: elastic shell.** Parameters `k_e` (the same for all particles) and `c_form`.

```
sig_ext^(i) = -k_e (u_surf^(i) - u_free),   u_free = c_form - c_ref
            = -K_eff (c_avg^(i) - c_form),   K_eff = k_e K / (k_e + K),   K = Lambda (1+3 nu_bar)
```

The shell is stress-free when the particle is at the filling `c_form` at which the shell formed (by default the initial filling). Substituting `u_surf = sig_ext/K + c_avg - c_ref` shows that the particle and the shell act as two springs in series. This holds for single- and mixed-phase particles alike, because `u_surf` depends only on `c_avg`. `K_eff` is stored as a derived parameter: it tends to `K` for a rigid medium (`k_e = Inf`) and to 0 when `k_e = 0` or `Lambda = 0`. The stress is integrated as `d sig_ext/dt = -K_eff (dc_avg/dt)/c_m`, starting from `-K_eff (c_0 - c_form)`, which keeps the relation exact (`elastic_stress`, `solve_stress_rate`).

Because the stress depends only on `c_avg`, the elastic shell is conservative: it adds `+K_eff (c_avg - c_form)` to every particle's potential on both branches, tilting the voltage curve by `K_eff R_g T/F` per unit filling, but it opens no hysteresis gap. It does change the ensemble: the slope of a particle's potential becomes `mu_chem' + K_eff` (single phase) or `K_eff - 2 A_int` (phase separated), so a shell with `K_eff > -min mu_chem'` (2 for the paper's `mu_chem`) stops particles from lithiating one after another. The binodal and spinodal points are unchanged.

**`'none'`:** `sig_ext = 0`.

To add another medium, extend `external_stress_law` and the matching branches of `solve_stress_rate`, `dae_residual` and `dae_jacobian`.

## Time integration

The ensemble model is an index-one DAE: `Phi` and the reaction currents `j_tr` are algebraic variables fixed by current conservation and Butler-Volmer kinetics, and the external-stress law is implicit in `d sig_ext/dt`. Each of these implicit parts has a unique solution for any state:

- the total current increases strictly with `Phi`, so `Phi` follows from a bracketed root solve (`solve_potential`);
- `j_tr` then follows from Butler-Volmer;
- the stress-law residual increases strictly with `d sig_ext/dt`, so the stress rate follows from bisection (`solve_stress_rate`).

By default (`solver = 'ode15s'`), the code therefore integrates only `y = [c_avg; sig_ext]` with `ode15s`, and the right-hand side (`reduced_rhs`) solves for `Phi`, `j_tr` and `d sig_ext/dt` at every evaluation. This is the same model, not an approximation: a test checks that the reduced right-hand side satisfies every equation of the DAE. It is more robust than the fully implicit form, because the steep switch of the plastic-shell stress is handled by bisection rather than by the integrator's Newton iteration. Phase changes are located with event functions; after each event the phase state is updated and integration resumes from the same state, with no consistent initial conditions needed. `Phi` and `j_tr` are rebuilt at the output times, so the solution has the same fields as before.

The fully implicit DAE in `(c_avg, j_tr, sig_ext, Phi)` is still available with `solver = 'ode15i'` (`solve_block`, `dae_residual`, `dae_jacobian`, `consistent_initial_conditions`). It can stop at current reversals with the surface-velocity stress model when many particles are present.

See [`docs/solver_formulation.md`](docs/solver_formulation.md) for the full description, the numerical settings and a suggested wording for the manuscript.

## Using your own chemical potential

The homogeneous chemical potential `mu_chem(c)` sets the material: its shape decides whether, and between which concentrations, a particle phase-separates. Everything else (critical points, the critical coupling `Lambda_c`, the single-particle curves and the ensemble) follows from it automatically.

### Where to change it

`mu_chem` is defined in one file, `src/thermodynamics/chemical_free_energy.m`, together with its integral and derivative:

```matlab
function [f, mu, dmu] = chemical_free_energy(c)
% c: normalised concentration in (0, 1); all outputs in units of R_g*T
f   = ...;      % g_chem(c), the homogeneous free energy
mu  = ...;      % mu_chem(c)  = d g_chem / dc
dmu = ...;      % mu_chem'(c) = d mu_chem / dc
end
```

`f_chem`, `mu_chem` and `mu_chem_prime` are thin wrappers around this file, so the rest of the code always sees the same function. The shipped version is the regular solution of the paper, `mu_chem = log(c/(1-c)) + 3(1-2c)`. To try another regular solution, change `Omega`; for a different form, replace the three expressions.

### Requirements

- **Consistency.** `mu` must be the derivative of `f`, and `dmu` the derivative of `mu`. `f` is used in the Gibbs minimisation, `mu` in the potentials and `dmu` in the stability analysis, so an inconsistent set gives wrong binodal or spinodal points. `test_mu_is_derivative_of_f` and `test_mu_prime_is_derivative_of_mu` check this by finite differences.
- **Vectorised and finite on (0, 1).** All three are evaluated element-wise on vectors of concentrations strictly between 0 and 1.
- **One miscibility gap.** Phase separation needs `mu_chem` to be non-monotonic. The model has three phase states (single phase, lithium-rich shell, lithium-poor shell), so it supports one two-phase region: `mu_chem'(c) + 2*A_int` must change sign once on each side of the global minimum of `mu_chem'`. Materials with several plateaus (e.g. graphite) need more phase states and are not covered.
- **No symmetry is assumed.** `mu_chem` need not be symmetric about `c = 1/2`.

### What follows automatically

- **Critical coupling.** Phase separation is possible only while `2*A_int < -min mu_chem'`, so `params.Lambda_c = -min(mu_chem')*3*(1+nu_bar)/(2*(1+3*nu_bar))` is computed from your `mu_chem` (`min_mu_chem_prime`). It is 1.8 for the paper's `mu_chem` and `nu = 0.25`, and 0 if `mu_chem` is monotonic. Above `Lambda_c` particles stay single phase as a solid solution (see [Binodal points](#binodal-points-analytical-vs-grid-search)).
- **Spinodal points** `c_sp`: the roots of `mu_chem'(c) + 2*A_int = 0` on either side of the global minimum of `mu_chem'`.
- **Binodal points** `c_b`: from the Gibbs-energy minimisation, evaluated at `c_avg = mean(c_sp)`, which always lies between the binodal points, and started from a guess built from `c_sp`. Override with `binodal_c_avg` and `binodal_initial_guess` if needed.

### Checking a new `mu_chem`

```matlab
setup_paths();
p = compute_critical_points(default_parameters());       % Lambda = 0
[p.Lambda_c, p.c_sp, p.c_b]
[c_b, info] = find_binodal_points(p, 'Method', 'both', 'Verbose', true);   % analytical vs grid search
run('examples/example_06_critical_points_vs_Lambda.m')   % critical points against Lambda
run('examples/example_05_ocv_hysteresis.m')              % single-particle curves
```

With `'Method', 'both'` the analytical binodal points are compared with a brute-force grid search; a warning means they disagree. Example 6 sweeps `Lambda` up to 1.6, so set its range below your `Lambda_c`.

Choose `soc_waypoints(1)` outside the spinodal range: `run_simulation` stops with `MechoCPM:unstableStart` if the cycle starts single phase between `c_s1` and `c_s2`, because the phase-switching rule is not defined there. Starting in the metastable range between a binodal and a spinodal point is fine.

### Tests that are specific to the paper's `mu_chem`

Most tests check relations that hold for any `mu_chem`. The following use the paper's `mu_chem` or its numbers, and are expected to fail after a change; adapt or skip them:

| Test | Uses |
|---|---|
| `test_free_energy_matches_paper`, `test_binodal_matches_closed_form` | The closed form of the regular solution with `Omega = 3` |
| `test_critical_coupling_Lambda_c`, `test_derived_mechanical_parameters` | `Lambda_c = 1.8`, symmetry about 1/2 |
| `test_regression_against_original_code` | Results saved from the original code |
| `test_ocv_without_phase_separation`, `test_cycle_without_phase_separation`, the `Lambda = 2` cases in `test_plotting` | `Lambda = 2, 3` being above `Lambda_c` |
| `test_unstable_initial_state_errors` | The spinodal points at `Lambda = 0.99` |
| `test_examples` | Example parameter values (`Lambda` up to 2) |

## Repository layout

```
MechoCPM/
|-- setup_paths.m                   adds src/ to the MATLAB path
|-- src/
|   |-- config/                     default_parameters, set_parameters, finalise_parameters
|   |-- solver/                     run_simulation, solve_block_ode15s, reduced_rhs, solve_potential,
|   |                               phase_change_events, update_phase_state, ramp_current,
|   |                               build_particle_population, build_current_protocol,
|   |                               solve_block, dae_residual, dae_jacobian, consistent_initial_conditions (ode15i)
|   |-- thermodynamics/             chemical_free_energy (+ f_chem, mu_chem, mu_chem_prime),
|   |                               mu_hys, find_binodal_points, find_spinodal_points, min_mu_chem_prime,
|   |                               compute_critical_points
|   |-- kinetics/                   exchange_current
|   |-- mechanics/                  external_stress_law, sig_ideal, surface_velocity, solve_stress_rate,
|   |                               elastic_stress, medium_label
|   |-- analysis/                   run_parameter_sweep, compute_ocv_hysteresis, voltage_summary
|   `-- plotting/                   plot_voltage_soc, plot_voltage_time, plot_potential, plot_external_stress,
|                                   plot_phase_fraction, plot_filling_ratio, plot_ocv_hysteresis, plot_gibbs_minimisers,
|                                   plot_results_comparison, new_figure, make_fixed_figure
|-- docs/solver_formulation.md      how the model is integrated in time
|-- examples/                       example_01 ... example_08, run_all_examples
|-- tests/                          test_thermodynamics, test_solver, test_examples, test_plotting, data/original_code_reference.mat
`-- .github/workflows/tests.yml     runs the tests on every push
```

## Examples

| Script | What it shows |
|---|---|
| `example_01_single_cycle` | One cycle with mechanics and external stress; all standard plots |
| `example_02_Lambda_sweep` | Effect of the mechanical coupling `Lambda` |
| `example_03_sigma_y_sweep` | Effect of the yield stress `sigma_y` |
| `example_04_crate_and_particles` | C-rate sweep, number of size groups, per-segment C-rates |
| `example_05_ocv_hysteresis` | Single-particle solution branches and cycling hysteresis with metastable and unstable states, no time integration (Fig. 3) |
| `example_06_critical_points_vs_Lambda` | Binodal (analytical and grid search) and spinodal points against `Lambda` (Fig. 2) |
| `example_07_external_stress_comparison` | Ensemble with and without external stress (Figs. 4-7) |
| `example_08_elastic_shell` | Elastic external medium: stiffness sweep, stress against filling, quasi-static curves |

## Output

`run_simulation` returns a struct with the time `t`, state of charge `SOC`, electrode potential `Phi`, applied current `I` and, per size group, `c_av`, `theta`, `jtr`, `sig_ext`, `mu` and phase state `s`, together with the `params` and `protocol` used, the number of phase-change events `n_events` and the `runtime`. With `verbose` on (the default) it also logs the electrode potential at the start and end of each block and its range over the run.

For the electrode potential:

| Function | Output |
|---|---|
| `plot_voltage_soc(sols, labels)` | `Phi` [mV] against SoC |
| `plot_voltage_time(sols, labels)` | `Phi` [mV] against time, with the protocol blocks marked |
| `voltage_summary(sols, labels, soc_points)` | Prints, and returns, `Phi` on the lithiation and delithiation branches at the given SoC (default 0.25, 0.5, 0.75), the gap between them and the range of `Phi` |
| `compute_ocv_hysteresis` | Single-particle branches (single phase with stable/metastable/unstable parts, mixed-phase line) and the lithiation and delithiation paths, as `mu` and `Phi` |
| `plot_gibbs_minimisers` | Fig. 3a,b: single- and mixed-phase branches |
| `plot_ocv_hysteresis(ocv)`, `plot_ocv_hysteresis(ocv, [], 'Phi')` | Fig. 3c,d: cycling paths as `mu`, or as the equilibrium potential `-mu*R_g*T/F` [mV] |

## Publishing this repository on GitHub

1. **Install Git** (not currently installed on this machine): download "Git for Windows" from <https://git-scm.com/download/win>, then open a new PowerShell window and check `git --version`. Tell Git who you are, once:
   ```bash
   git config --global user.name "Your Name"
   git config --global user.email "you@example.com"
   ```
2. **Create the repository on GitHub.** Sign in at <https://github.com>, click **New repository**, choose a name (for example `MechoCPM`), set it to **Public** (or Private for now), and leave "Add a README", ".gitignore" and "license" **unticked**: they already exist locally. Click **Create repository** and copy its URL, e.g. `https://github.com/<your-username>/MechoCPM.git`.
3. **Choose a licence** and add it as a file named `LICENSE` in this folder, and replace the `TODO` in the section below. MIT or BSD-3-Clause are common for academic code. Without a licence, others have no legal right to reuse the code.
4. **Push the code** from PowerShell, inside the folder that contains this README:
   ```bash
   cd path\to\MechoCPM
   git init
   git add .
   git status
   git commit -m "Initial release: code accompanying the JES paper"
   git branch -M main
   git remote add origin https://github.com/<your-username>/MechoCPM.git
   git push -u origin main
   ```
   `git status` before committing should list the source files, `tests/data/original_code_reference.mat` and no stray `.mat` or `.asv` files (`.gitignore` excludes them). GitHub may ask you to sign in in a browser the first time.
5. **Check the tests on GitHub.** Open the repository's **Actions** tab: the "MATLAB tests" workflow runs on every push and should turn green after a few minutes. If it asks for a MATLAB licence, check that the repository is public (free for public repositories).
6. **Finish the repository page.** Add a short description and topics in the repository's settings (About box), and replace the `doi` placeholder in `CITATION.cff` once the paper is published.
7. **Optional: get a citable DOI.** Link the repository to <https://zenodo.org>, then create a GitHub *Release* (for example `v1.0.0`). Zenodo archives it and gives it a DOI that you can cite in the paper's "Code availability" statement.

If you prefer not to use the command line, **GitHub Desktop** (<https://desktop.github.com>) does steps 1, 2 and 4 with buttons: *File > Add local repository*, then *Publish repository*.

## Citation

If you use this code, please cite the paper (see `CITATION.cff`).

## License
MIT License
