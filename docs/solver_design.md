# Solver design

## 1. Objective and the Levenberg–Marquardt step

Minimise `F(x) = ½ |r(x)|²` over the point coordinates `x ∈ ℝⁿ`, with
`r : ℝⁿ → ℝᵐ` the concatenated constraint residuals and `J = ∂r/∂x` (m × n).

**Gauss–Newton model.** Linearise `r(x + δ) ≈ r + Jδ`. Then
`F(x + δ) ≈ ½|r|² + δᵀJᵀr + ½ δᵀJᵀJ δ`. Setting the gradient to zero gives the
Gauss–Newton step `JᵀJ δ = −Jᵀr`. Two things can go wrong: `JᵀJ` is singular
for under-constrained systems (many valid `δ`, or none that is stable), and the
linearisation may be poor far from a solution so the step overshoots.

**Damping.** Levenberg–Marquardt minimises the model plus a trust penalty,
`½|r + Jδ|² + ½ λ |δ|²`, whose stationarity condition is

```
(JᵀJ + λ I) δ = −Jᵀr .
```

For `λ > 0` the matrix is symmetric positive definite whatever the rank of `J`
(`vᵀ(JᵀJ + λI)v ≥ λ|v|² > 0`), so a unique, finite step always exists. Limits:
`λ → 0` recovers Gauss–Newton (fast, quadratic for zero-residual problems);
`λ → ∞` gives `δ ≈ −Jᵀr/λ`, a short steepest-descent step that always decreases
`F` for large enough `λ`. Adapting `λ` interpolates between the two.

The system is solved with `Eigen::LDLT` (a Cholesky-type factorisation of an
SPD matrix). **No matrix inverse is formed.** Forming `JᵀJ` squares the
condition number; with `λ` added this is acceptable at the sizes this library
targets, and the alternative (QR on the augmented matrix `[J; √λ I]`) is listed
as a future improvement.

## 2. Algorithm

```
validate options and problem            → InvalidInput on failure
x ← encode(initial points); r ← r(x)    → InvalidInput if invalid/non-finite
λ ← initialDamping
loop:
    if max|r| ≤ residualTolerance       → Converged
    if iterations ≥ maxIterations       → MaxIterations
    J ← finite-difference Jacobian(x)    → NumericalFailure if it throws / non-finite
    g ← Jᵀ r
    if max|g| ≤ gradientTolerance       → Stagnated   (stationary, residual still large)
    repeat:
        if λ > maxDamping               → Stagnated or NumericalFailure
        solve (JᵀJ + λI) δ = −g with LDLT
            failure or non-finite δ     → λ ← λ·increase; reject; continue
        x_trial ← x + δ; evaluate r(x_trial)
            throws / non-finite         → treat as rejected
        ρ ← (F(x) − F(x_trial)) / (½ δᵀ(λδ − g))      ("gain ratio")
        if F decreased and ρ ≥ minGainRatio:
            accept: x ← x_trial, r ← r_trial, λ ← max(λ·decrease, minDamping)
            if ‖δ‖ ≤ stepAbs + stepRel‖x‖ and max|r| > tol  → Stagnated
        else:
            reject: λ ← λ·increase
```

Every exit goes through one `finish` routine that fills the final state, the
per-constraint residuals and the Jacobian analysis, so *every* status carries
diagnostics.

### Why the acceptance test is justified

For the damped model the predicted decrease is
`m(0) − m(δ) = −gᵀδ − ½δᵀJᵀJδ`. Using the normal equation `JᵀJδ = −g − λδ`
this equals `½ δᵀ(λδ − g)`, which is strictly positive whenever `δ ≠ 0`
(`−gᵀδ = δᵀ(JᵀJ + λI)δ > 0`). The ratio `ρ = actual / predicted` therefore
measures how trustworthy the linear model was. Requiring a *positive actual
decrease* guarantees `F` is monotonically non-increasing — the property the tests
check for the cost history; the small `minGainRatio` (1e-4) rejects steps whose
improvement is dominated by rounding in the predicted decrease. A rejected step
increases `λ`, shrinking the next step and rotating it toward steepest descent,
for which a decrease is guaranteed once `λ` is large enough (up to floating-point
resolution of `F`, which is why `maxDamping` exists as a termination criterion).

### Damping strategy

* Accepted step: `λ ← max(0.3 λ, minDamping)` — trust the Gauss–Newton model more.
* Rejected step: `λ ← 4 λ` — tighten the trust region. Failed linear solves
  and degenerate trial states are treated like rejections.
* Defaults `λ₀ = 1e-3`; `minDamping = 1e-14`, `maxDamping = 1e12`.

The identity damping (instead of Marquardt's `diag(JᵀJ)` scaling) is what the
specification asks for. It is not scale invariant; the residuals here are
lengths or unit-direction quantities and the Jacobian columns are O(1), so the
unscaled form behaves well. See `limitations.md` for the consequence for
badly scaled custom constraints.

## 3. Finite-difference Jacobian

Column `j` perturbs coordinate `xⱼ` only. With step
`hⱼ = h·max(1, |xⱼ|)` (`h = finiteDifferenceStep`, default `1e-6`):

```
central :  ∂r/∂xⱼ ≈ ( r(x + hⱼeⱼ) − r(x − hⱼeⱼ) ) / (2 hⱼ)      error O(h²·r''')
forward :  ∂r/∂xⱼ ≈ ( r(x + hⱼeⱼ) − r(x) ) / hⱼ                  error O(h·r'')
```

* **Why scale by `max(1,|xⱼ|)`:** a perturbation must be large compared with
  the spacing of doubles near `xⱼ` (`ε|xⱼ|`, ε ≈ 2.2e-16) to avoid cancellation,
  yet small compared with the length scale over which the residual curves. The
  relative scale handles large coordinates; the floor `1` avoids a vanishing step
  near the origin.
* **Why 1e-6 and central by default:** the rounding error of a difference
  quotient is ≈ `ε|r|/h` and the central truncation error is ≈ `h²|r'''|/6`.
  Balancing them puts the optimum near `h ≈ ε^(1/3) ≈ 6e-6`, giving ≈1e-10 – 1e-11
  accuracy; the observed agreement with analytic gradients in
  `constraint_tests.cpp` is better than 1e-8 absolute for coordinates of order
  1–15. Forward differences need `h ≈ √ε ≈ 1.5e-8` for the same balance, but
  only reach ≈1e-8 – 1e-7 accuracy.
* The step actually used is the representable difference `(x+h)−x`, removing
  the rounding bias of the perturbation itself.
* The finite-difference Jacobian is built column by column with a full
  residual evaluation each (2n evaluations for central); the Jacobian is dense.
  This is simple and adequate for small systems; analytic or sparse
  Jacobians are a listed extension.

`finiteDifferenceJacobian` is public and verified against analytic derivatives
for the distance, plane, coincidence, angle and parallel constraints.

## 4. Convergence and termination

| Status | Trigger | Meaning |
|---|---|---|
| `Converged` | `max|rᵢ| ≤ residualTolerance` | all constraints satisfied to tolerance. **Only** this status means success. |
| `MaxIterations` | accepted iterations ≥ `maxIterations` | not converged yet; the state is the best iterate so far |
| `Stagnated` | `max|Jᵀr| ≤ gradientTolerance`, or accepted step `≤ stepAbs + stepRel‖x‖`, or no damped step reduces `F` (`λ > maxDamping`) — in each case with the residual still above tolerance | a stationary point of `F` with nonzero residual: inconsistent constraints or a local minimum (including symmetric saddle points) |
| `InvalidInput` | bad options, empty problem, non-finite initial coordinates, a constraint that cannot be evaluated at the initial state (e.g. zero-length direction) | nothing was solved |
| `NumericalFailure` | the Jacobian cannot be evaluated or is non-finite mid-solve; LDLT fails for every damping value | the iteration broke down |

Rules:

* `residualTolerance` is **absolute** and applied to the max-norm of the
  residuals (units: length for distance/point/plane constraints, dimensionless
  for direction constraints). A single relative tolerance would be ambiguous for a
  vector mixing both kinds of quantities, and a tolerance relative to the
  *initial* residual would declare success on any large improvement.
* The step test mixes an absolute and a relative part (`stepAbs + stepRel‖x‖`),
  because coordinates may be large.
* The gradient test uses an absolute value; `∥Jᵀr∥` scales with both the
  residual and the Jacobian.
* Reaching the iteration loop's end is never reported as success.

## 5. Options

| Option | Default | Role |
|---|---|---|
| `maxIterations` | 100 | cap on accepted iterations |
| `residualTolerance` | 1e-9 | success criterion on `max|r|` |
| `stepAbsTolerance`, `stepRelTolerance` | 1e-12, 1e-12 | stagnation by tiny accepted step |
| `gradientTolerance` | 1e-10 | stagnation at a stationary point |
| `initialDamping` | 1e-3 | initial `λ` |
| `dampingDecrease`, `dampingIncrease` | 0.3, 4 | `λ` multipliers on accept / reject |
| `minDamping`, `maxDamping` | 1e-14, 1e12 | clamp and give-up thresholds |
| `minGainRatio` | 1e-4 | minimum `ρ` to accept |
| `finiteDifferenceStep` | 1e-6 | relative step `h` |
| `jacobianScheme` | `Central` | `Forward` or `Central` |
| `rankTolerance` | 1e-8 | relative singular-value cutoff for the diagnostics |

`SolverOptions::validate()` rejects non-finite or out-of-range values, and the
solver turns that into `InvalidInput` rather than throwing.

## 6. Diagnostics

Every `SolverResult` contains: status and message; final points (the best
iterate, which is only a solution when `Converged`); accepted/rejected step
counts and residual-evaluation count; initial/final residual norms, final
max-norm and final gradient norm; final damping; the cost history (strictly
decreasing); per-constraint residual vectors with their norms; and a
`SystemAnalysis` with the Jacobian's numerical rank, local DOF, redundant
residual count, rank threshold and singular values. `formatReport` renders all
of this as text.

**Convergence is not classification.** A converged system can still have free
degrees of freedom (`Converged` with local DOF = 2 in the sphere example),
redundant constraints can be consistent (`Converged` with redundant rows > 0),
and inconsistent ones end in `Stagnated`. The two pieces of information are
reported separately; see `numerical_robustness.md` §4 for what the DOF estimate
does and does not establish.
