# MechoCPM

**Chemo-mechanical model of phase-changing electrode particles**

MATLAB code accompanying *Coupled chemo-mechanical framework for phase-changing electrodes* (J. Petrassem de Sousa, W. Clarke, A. Galvis, S. Sahu and J. M. Foster, Journal of The Electrochemical Society, 2026).

The code simulates galvanostatic cycling of an ensemble of electrode particles with a distribution of sizes. Each particle is either single-phase or phase-separated into a core and a shell. Transitions are triggered when the particle's average concentration `c_avg` crosses a spinodal point (single phase -> core-shell) or a binodal point (core-shell -> single phase). Mechanics enters in two independent ways, following the scaling of the paper. The mechanical coupling `Lambda`, the ratio of the elastic misfit energy of a fully lithiated particle to the thermal energy, sets the internal incompatibility stress between core and shell: it shifts the spinodal and binodal points and tilts the mixed-phase branch of the potential. An external radial stress `sigma_ext`, which relaxes towards a plastic yield stress `sigma_y`, translates the potential rigidly. Stresses are measured in units of `R_g*T/V_c` (about 0.3 GPa), so the external stress enters the potential with unit coefficient.

Symbols in the code follow the paper (see [Notation](#notation)); equation numbers in the source comments refer to the paper.

## Contents

- [Quick start](#quick-start)
- [Testing that everything works](#testing-that-everything-works)
- [Choosing parameters](#choosing-parameters)
- [Binodal points: analytical vs grid search](#binodal-points-analytical-vs-grid-search)
- [Notation](#notation)
- [Model summary](#model-summary)
- [Using your own free energy](#using-your-own-free-energy)
- [Repository layout](#repository-layout)
- [Examples](#examples)
- [Output](#output)
- [Publishing this repository on GitHub](#publishing-this-repository-on-github)

## Quick start

Requirements: MATLAB R2019b or newer. The Optimization Toolbox (`fsolve`) is needed for the default analytical binodal solver; `ode15i` and `fzero` are in base MATLAB, and the grid-search binodal method needs no toolbox.

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

All 33 tests should pass (`test_examples` runs examples 1-6 without figures on screen, so allow a minute or two). They check, among other things, that

- `mu_chem`, `mu_chem_prime` and `f_chem` are mutually consistent and equal the paper's `mu_chem`;
- the spinodal points satisfy the coherent-spinodal condition `mu_chem'(c) + 2*A_int = 0`, exist only for `Lambda < Lambda_c = 1.8`, and the binodal points satisfy the common-tangent construction at `Lambda = 0`;
- the analytical and grid-search binodal methods agree, and the binodal points do not depend on `c_avg`;
- `mu_hys` equals the paper's `mu_single` and `mu_sep`, the core and shell potentials are equal at equilibrium, and the external stress shifts the potential by exactly `-sig_ext` for any `Lambda`;
- old `gamma` inputs are rejected with the conversion to `Lambda`, and examples 1-6 run and produce figures and logs that never mention `gamma`;
- a hysteretic external stress of `+-sigma_y` widens the simulated voltage gap by exactly `2*sigma_y*R_g*T/F`, and the quasi-static equilibrium potential is `-mu*R_g*T/F`;
- the DAE initial conditions are consistent and a short cycle conserves charge;
- **`test_regression_against_original_code`**: the electrode potential, binodal and spinodal points reproduce results saved from the original, pre-reorganisation code (`tests/data/original_code_reference.mat`) for four parameter sets, to within solver tolerance (RMS difference below 10 microvolts). The reference was made with the previous scaling and is converted with `Lambda = gamma*alpha`, `sigma_y = gamma*sigma_y_old`, which checks that the change of scaling leaves the model unchanged.

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
| `example_03_sigma_y_sweep` | seconds | Electrode potential against SoC and time for `sigma_y` = 0, 0.2, 0.4, 0.8: larger `sigma_y` widens the gap between lithiation and delithiation (about 20 mV at `sigma_y` = 0 to about 60 mV at 0.8 near SoC = 0.5). Prints the simulated increase of the gap next to the prediction `2*sigma_y*R_g*T/F` (about 52 mV per unit); the two agree |
| `example_04_crate_and_particles` | seconds | Electrode potential against SoC and time, with printed summaries, for a C-rate sweep, a number-of-groups sweep and per-segment C-rates |
| `example_05_ocv_hysteresis` | 1 s | Paper Fig. 3a (flat plateau at `Lambda = 0`), Fig. 3b (sloped mixed-phase line, `Lambda = 0.99`) and Fig. 3d (curves shifted by `sig_ext = -0.6, 0, 0.6`), each with the matching equilibrium electrode potential `-mu*R_g*T/F` [mV]. Prints `c_b = [0.0707 0.9293]` for `Lambda = 0` and `c_b = [0.1570 0.8430]`, `c_sp = [0.2857 0.7143]` for `Lambda = 0.99`, and the potentials at c_avg = 0.5, which shift by about +-15.5 mV for `sig_ext = -+0.6` |
| `example_06_critical_points_vs_Lambda` | seconds | Paper Fig. 2: spinodal and binodal points against `Lambda`, with `Lambda_c = 1.8` marked; prints the largest analytical vs grid-search difference (below 1e-3) |
| `example_07_external_stress_comparison` | 1-2 min | Paper Figs. 4/6 (filling ratio of each particle group), Fig. 7 (electrode potential against SoC with and without external stress, `Lambda = 1.32`, `sigma_y = 0.8`) and the electrode potential against time. Prints both voltage summaries and the gap increase, about 41 mV (`2*0.8*R_g*T/F`) |

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
| `sigma_y` | 0 | sigma_y^* | Scaled yield stress of the external shell, `V_c*sigma_y/(R_g*T)` |
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
| `k` | 0.01 1/s | Relaxation rate of the external stress |
| `u0` | 1e-5 A | Current scale in sigma_ideal |
| `I0` | 1 A/m^2 | Exchange current density |
| `c_m` | 2.7e5 mol/m^3 | Maximum concentration |
| `F`, `R_g`, `T` | 96500, 8.314, 300 | Faraday constant, gas constant, temperature |
| `compute_critical_points` | `true` | If `false`, `c_b` and `c_sp` are used as given |
| `binodal_method` | `'analytical'` | `'analytical'`, `'gridsearch'` or `'both'` (see below) |
| `binodal_grid_points`, `binodal_tol` | 2000, 1e-3 | Grid resolution and agreement tolerance for the grid search |
| `binodal_initial_guess`, `binodal_c_avg` | `[0.1 0.9]`, 0.51 | `fsolve` starting point and c_avg used in the minimisation |
| `RelTol`, `AbsTol` | 1e-6, 1e-8 | `ode15i` tolerances |
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

The spinodal points `c_sp = [c_s1, c_s2]` come from the coherent-spinodal condition `mu_chem'(c) + 2*A_int = 0` in `find_spinodal_points`. For `Lambda >= Lambda_c` (1.8 for the default `mu_chem` and `nu = 0.25`) no spinodal points exist: phase separation is suppressed and `run_simulation` stops with an explanatory error. The external stress does not affect either set of critical points.

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
| k, u0, I0, c_m | `k`, `u0`, `I0`, `c_m` | As in the ensemble parameter table |

The phase state `s` of a particle is `0` (single phase), `+1` (lithium-rich shell, core `c_b1`, shell `c_b2`) or `-1` (lithium-poor shell, core `c_b2`, shell `c_b1`).

## Model summary

For size group `i` (radius `R_i`) with average concentration `c_avg^(i)`, phase state `s_i`, external stress `sig_ext^(i)` and a common electrode potential `Phi`:

```
dc_avg^(i)/dt  = -j_tr^(i) A^(i) / (F V^(i))
eta^(i)        = (F Phi + mu^(i)) / (R_g T)
j_tr^(i)       = 2 I0 sqrt(c_surf (c_m - c_surf)) / c_m * sinh(eta/2)
I(t)           = sum_i N^(i) A^(i) j_tr^(i)
d sig_ext/dt   = -k (sig_ext - sig_ideal),
sig_ideal      = sigma_y tanh(I/u0)
```

`mu^(i)` is the surface chemo-mechanical potential from the single-particle model (`mu_hys`): in units of `R_g*T`,

```
mu_single = mu_chem(c_avg) - sig_ext
mu_sep    = mu_chem(c_s) - 2 A_int (c_avg - c_s) - sig_ext,   A_int = Lambda (1+3 nu_bar) / (3 (1+nu_bar))
```

for a single-phase and a core-shell particle (shell concentration `c_s`). `Lambda` tilts the mixed-phase branch; the external stress translates both branches by `-sig_ext`, so a hysteretic stress of `+-sigma_y` widens the voltage gap by `2*sigma_y*R_g*T/F`. The system is a fully implicit DAE solved with `ode15i`. Phase changes are located with event functions; after each event, consistent initial conditions are recomputed and integration resumes. At the start of each current block the external stress is reset to `sign(I)*sigma_y` (the shell is assumed to have yielded in the direction of the new current), and the current is ramped smoothly between blocks (`ramp_current`).

## Using your own free energy

The chemical free energy is defined in a single file, `src/thermodynamics/chemical_free_energy.m`, which returns `g_chem`, `mu_chem = g_chem'` and `mu_chem'`. `f_chem`, `mu_chem` and `mu_chem_prime` are thin wrappers around it, so the three can never become inconsistent. The shipped version is the regular solution used in the paper, `mu_chem = log(c/(1-c)) + 3(1-2c)`. Replace it with your own expression and run `runtests('tests')`: the tests check that the derivatives are consistent and that the critical points satisfy the common-tangent and stability conditions (the regression test against the original results will, of course, no longer apply).

## Repository layout

```
MechoCPM/
|-- setup_paths.m                   adds src/ to the MATLAB path
|-- src/
|   |-- config/                     default_parameters, set_parameters, finalise_parameters
|   |-- solver/                     run_simulation, solve_block, dae_residual,
|   |                               phase_change_events, consistent_initial_conditions,
|   |                               ramp_current, build_particle_population, build_current_protocol
|   |-- thermodynamics/             chemical_free_energy (+ f_chem, mu_chem, mu_chem_prime),
|   |                               mu_hys, find_binodal_points, find_spinodal_points,
|   |                               compute_critical_points
|   |-- kinetics/                   exchange_current
|   |-- mechanics/                  sig_ideal
|   |-- analysis/                   run_parameter_sweep, compute_ocv_hysteresis, voltage_summary
|   `-- plotting/                   plot_voltage_soc, plot_voltage_time, plot_potential, plot_external_stress,
|                                   plot_phase_fraction, plot_filling_ratio, plot_ocv_hysteresis, plot_ocv_voltage,
|                                   plot_results_comparison, new_figure, make_fixed_figure
|-- examples/                       example_01 ... example_07, run_all_examples
|-- tests/                          test_thermodynamics, test_solver, test_examples, data/original_code_reference.mat
`-- .github/workflows/tests.yml     runs the tests on every push
```

## Examples

| Script | What it shows |
|---|---|
| `example_01_single_cycle` | One cycle with mechanics and external stress; all standard plots |
| `example_02_Lambda_sweep` | Effect of the mechanical coupling `Lambda` |
| `example_03_sigma_y_sweep` | Effect of the yield stress `sigma_y` |
| `example_04_crate_and_particles` | C-rate sweep, number of size groups, per-segment C-rates |
| `example_05_ocv_hysteresis` | Quasi-static single-particle potential, no time integration (Fig. 3) |
| `example_06_critical_points_vs_Lambda` | Binodal (analytical and grid search) and spinodal points against `Lambda` (Fig. 2) |
| `example_07_external_stress_comparison` | Ensemble with and without external stress (Figs. 4-7) |

## Output

`run_simulation` returns a struct with the time `t`, state of charge `SOC`, electrode potential `Phi`, applied current `I` and, per size group, `c_av`, `theta`, `jtr`, `sig_ext`, `mu` and phase state `s`, together with the `params` and `protocol` used, the number of phase-change events `n_events` and the `runtime`. With `verbose` on (the default) it also logs the electrode potential at the start and end of each block and its range over the run.

For the electrode potential:

| Function | Output |
|---|---|
| `plot_voltage_soc(sols, labels)` | `Phi` [mV] against SoC |
| `plot_voltage_time(sols, labels)` | `Phi` [mV] against time, with the protocol blocks marked |
| `voltage_summary(sols, labels, soc_points)` | Prints, and returns, `Phi` on the lithiation and delithiation branches at the given SoC (default 0.25, 0.5, 0.75), the gap between them and the range of `Phi` |
| `compute_ocv_hysteresis` + `plot_ocv_voltage` | Quasi-static equilibrium potential `-mu*R_g*T/F` of a single particle |

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
