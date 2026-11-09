# Solver formulation: reduced ODE with `ode15s`

This note explains how MechoCPM integrates the ensemble model in time, why the default changed from a fully implicit DAE solved with `ode15i` to a reduced ODE solved with `ode15s`, and why the two describe the same model.

## 1. The model as a DAE

For `n` particle-size groups, the ensemble model of the paper ("Ensemble of interacting particles") has, for each group `i`,

```
dc_i/dt   = -j_i A_i / (F V_i)                                   (lithium conservation)
0         = BV(c_i, s_i, Phi, sig_i) - j_i                       (Butler-Volmer)
dsig_i/dt = -k (sig_i - sig_ideal(u_dot_i))                      (external stress, plastic shell)
u_dot_i   = (dsig_i/dt) / (Lambda (1+3 nu_bar)) + (dc_i/dt) / c_m
sig_ideal = -sigma_y tanh(u_dot / u0)
```

(for an elastic shell the stress row is instead `dsig_i/dt = -K_eff (dc_i/dt)/c_m`, see the README), and, for the whole electrode,

```
0 = I(t) - sum_i N_i A_i j_i                                     (current conservation)
```

Here `c_i` is the average concentration, `j_i` the reaction current density, `sig_i` the scaled external stress, `s_i` the phase state and `Phi` the electrode potential. `BV` is the Butler-Volmer law with the surface chemo-mechanical potential `mu_hys` (see `exchange_current`).

Two features make this system implicit:

- `j_i` and `Phi` are **algebraic** variables. They have no equation of their own and are fixed by the Butler-Volmer and current-conservation constraints.
- The stress law is **implicit in its own rate**: `dsig_i/dt` appears inside the `tanh` through `u_dot_i`, so it cannot be written as `dsig_i/dt = f(...)` in closed form.

## 2. The previous approach: fully implicit DAE with `ode15i`

The original code wrote all `3n + 1` equations as `F(t, w, w') = 0` with `w = [c; j; sig; Phi]` (`dae_residual`) and integrated them with `ode15i`. At every step, `ode15i` solves the whole system for `w` and `w'` together with Newton's method.

This works well while the equations are smooth. It ran into trouble with the external-stress model driven by the surface velocity. In that model, `sig_ideal` switches between `+sigma_y` and `-sigma_y` over a surface velocity of order `u0`, which is much smaller than the swelling rate. Wherever a particle's surface starts or stops moving, for example at every current reversal, the sensitivity of the stress equation to `dsig/dt` jumps by about four orders of magnitude within a fraction of a second. Newton's linearised steps overshoot there, the iteration fails to converge, and `ode15i` keeps reducing the step until it stops ("Unable to meet integration tolerances…"). An exact Jacobian (`dae_jacobian`) got the integration through lithiation, but not through the reversal with 10-20 particles at C/100.

## 3. The reduced ODE

The state passed to the integrator is only

```
y = [c_1 ... c_n, sig_1 ... sig_n]
```

and the right-hand side `reduced_rhs` computes `dy/dt` from `y` by solving the implicit parts itself:

1. **Electrode potential.** Solve current conservation, `sum_i N_i A_i BV_i(Phi) = I(t)`, for `Phi` (`solve_potential`). The total current increases strictly with `Phi`, so the root is unique. It is bracketed by stepping outwards from the previous `Phi` (a warm start) and refined with `fzero` to 1e-14 V.
2. **Reaction currents and concentration rates.** `j_i = BV(c_i, s_i, Phi, sig_i)`, then `dc_i/dt = -j_i A_i / (F V_i)`.
3. **Stress rates.** Solve the stress law for `dsig_i/dt` (`solve_stress_rate`). The residual `g(x) = x + k (sig_i + sigma_y tanh(u_dot_i(x)/u0))` increases strictly with `x` and changes sign on `[-k(sig_i + sigma_y), -k(sig_i - sigma_y)]`. The root is therefore unique, and bisection on that interval finds it to machine precision, however steep the `tanh` is. For the elastic shell the rate is explicit, `dsig_i/dt = -K_eff (dc_i/dt)/c_m`, and no root solve is needed.

`ode15s` then integrates `dy/dt = f(t, y)` as an ordinary (stiff) ODE.

### Why this is the same model

Each implicit part has a unique solution for any state `(c, sig)`, so the rates are a well-defined function of the state, even without a closed form. This is the same condition that makes the DAE index 1. `ode15i` also relies on it; it just finds the solution with Newton's method inside each step rather than with dedicated root finders. Nothing is approximated or eliminated beyond solving the same equations.

This is checked directly by `test_reduced_rhs_satisfies_dae_residual`: the rates and algebraic values returned by `reduced_rhs` make every row of `dae_residual` vanish.

### Why it is more robust

The steep stress switch is now handled by bisection, which cannot fail on a monotone function. The integrator only sees the resulting rates, which vary smoothly in time. The potential is found by a bracketed root solve rather than as part of a coupled Newton iteration.

## 4. Phase changes

Phase changes are located with the same event function as before (`phase_change_events`): a single-phase particle switches when its `c_avg` reaches a spinodal point, and a core-shell particle when it reaches a binodal point. At an event, `ode15s` stops, the phase state is updated (`update_phase_state`) and integration resumes from the same `y`.

Unlike `ode15i`, no consistent initial conditions are needed after an event or at the start of a block. `Phi` and `j` are not integrated, so they cannot become inconsistent with `(c, sig)`; they are recomputed at every evaluation.

## 5. Output

The solution has the same layout as before. At each output time, `solve_block_ode15s` rebuilds `j` and `Phi` from `(c, sig)` with `solve_potential`, so `sol.Phi`, `sol.jtr` and the other fields of `run_simulation` are unchanged.

## 6. Numerical settings

| Parameter | Default | Role |
|---|---|---|
| `solver` | `'ode15s'` | `'ode15s'` (reduced ODE) or `'ode15i'` (fully implicit DAE) |
| `RelTol`, `AbsTol` | 1e-6, 1e-8 | Integration tolerances (both solvers) |
| `max_step` | `[]` | Largest `ode15s` step [s]. Empty means 1/1000 of each block, so the output resolves the voltage curve to about 0.001 in SoC; otherwise `ode15s` takes very long steps on smooth stretches and the curve looks polygonal |
| `u0` | `[]` | Reference surface velocity in `sig_ideal` [1/s]. Empty means 5 % of the nominal swelling rate `min(c_rate)/3600` (stored as `u0_used`), which keeps the shell rigid-plastic (`tanh(20) = 1`) at any C-rate |

The inner solves are far more accurate than the integration tolerance (`Phi` to 1e-14 V, `dsig/dt` to machine precision), so they do not disturb `ode15s`'s step-size control. The Jacobian of `f` is formed by `ode15s` with finite differences.

## 7. Choosing a solver

```matlab
sol = run_simulation('Lambda', 1.32, 'sigma_y', 0.8);                     % ode15s (default)
sol = run_simulation('Lambda', 1.32, 'sigma_y', 0.8, 'solver', 'ode15i'); % fully implicit DAE
```

The `ode15i` path is kept for comparison and for cases where it is known to work. `test_ode15s_matches_ode15i` checks that the two give the same electrode potential (RMS difference below 0.1 mV), the same number of phase-change events and the same external stresses.

## 8. Validation

With `ode15s` as the default, the full test suite passes, including:

- the regression test against the original, pre-reorganisation code (electrode potential RMS difference below 10 microvolts, for the reference cases without external stress);
- all quick examples, including examples 1 and 3, which `ode15i` could not complete with the surface-velocity stress model;
- the external-stress checks: direction, continuity across blocks, the switching window `2 sigma_y / (Lambda (1+3 nu_bar))` and the gap `2 sigma_y R_g T / F`.

## 9. Cost

Each evaluation of `f` solves for `Phi` (a bracketed `fzero`, usually a few iterations thanks to the warm start) and runs a vectorised bisection for the stress rates. An evaluation is therefore more expensive than one residual evaluation of `ode15i`, but no Newton iterations on the full system are needed and there are no failed steps at the stress switches. The step cap `max_step` sets a minimum of 1000 steps per block.

## 10. Where things are

| File | Role |
|---|---|
| `src/solver/run_simulation.m` | Protocol loop; dispatches on `params.solver` |
| `src/solver/solve_block_ode15s.m` | One block with `ode15s`, events, output reconstruction |
| `src/solver/reduced_rhs.m` | Right-hand side of the reduced ODE |
| `src/solver/solve_potential.m` | `Phi` and `j` from the state (shared with `ode15i`) |
| `src/mechanics/solve_stress_rate.m` | `dsig/dt` from the implicit stress law |
| `src/solver/update_phase_state.m` | Phase-switching rule (shared) |
| `src/solver/solve_block.m`, `dae_residual.m`, `dae_jacobian.m`, `consistent_initial_conditions.m` | Fully implicit DAE path (`ode15i`) |
