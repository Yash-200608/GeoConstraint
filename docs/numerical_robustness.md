# Numerical robustness

## 1. Coordinate scale and the central tolerance policy

All geometric thresholds live in one struct, `Tolerances`
(`numerical/tolerances.hpp`); there is no other hidden constant in the geometry
layer, and the default instance is `constexpr` (no mutable global state). Every
geometry function that needs a threshold accepts a `Tolerances` argument.

| Field | Default | Used for |
|---|---|---|
| `zeroLength` | 1e-12 | a vector at most this long is "zero" and cannot define a direction |
| `length` | 1e-9 | absolute part of length comparisons (coincidence, containment) |
| `relative` | 1e-9 | relative part, multiplied by the magnitude of the operands |
| `angular` | 1e-9 | `sin θ` for parallelism, `cos θ` for perpendicularity (unit vectors) |

**Intended coordinate scale:** model units of roughly 1e-3 … 1e6 (millimetres to
kilometres). The absolute thresholds are meaningful in that range; the
`relative` term exists because doubles near 1e6 only resolve ≈1e-10, so a purely
absolute 1e-9 comparison would sit at the rounding level.

**Absolute + relative comparison.** `approxEqual(a,b)` accepts
`|a−b| ≤ length + relative·max(|a|,|b|)`. The absolute part protects values near
zero (where a purely relative test can never succeed); the relative part widens
the test for large values. The cost: at large magnitudes the test is loose in
absolute terms — two lines 1e-3 apart near x = 1e7 are reported `Coincident`
(tolerance ≈ 1e-2), while 1e-3 apart at the origin they are `Parallel`. Both
behaviours are asserted in `numerical_tests.cpp`. Callers who need different
semantics pass their own `Tolerances`.

## 2. Degenerate geometry is rejected, not repaired

| Situation | Behaviour |
|---|---|
| zero / non-finite direction passed to `normalized`, `Line3`, `Plane3`, `Matrix3::rotationAxisAngle`, `angleBetween`, `isParallel`, `isPerpendicular` | `DegenerateGeometryError` |
| non-finite origin / coordinates in a constructor | `InvalidArgumentError` |
| division of a vector by zero / non-finite scalar | `InvalidArgumentError` |
| plane from coincident or (nearly) collinear points: `|(b−a)×(c−a)| ≤ angular·|b−a||c−a|` | `DegenerateGeometryError` |
| a rotation matrix that is not proper | `InvalidArgumentError` from `RigidTransform3` |
| line–plane: parallel | `NoIntersection` / `LineInPlane` (never an arbitrary point) |
| line–line: parallel | `Parallel` / `Coincident`, distance still exact |
| non-finite results | `Invalid` relation / empty `optional`s |
| constraint evaluated with a zero direction | `DegenerateGeometryError`; inside `solve`, `InvalidInput` (at the start) or a rejected trial step (mid-solve) |

`Line3` and `Plane3` can only be created through validating factories, so a
live object is never degenerate.

**Near-degenerate is not degenerate.** A line tilted 1e-6 rad from a plane is a
legitimate unique crossing, but the crossing is ≈1e6 away from the line's origin
(`t = h/sin ε`); a tilt of 1e-12 rad is below `angular` and is classified as
parallel. The threshold is a *model decision*: very shallow crossings are
numerically valid but physically suspect. Lines tilted by `ε` have closest-point
parameters of size `1/ε`, so the *position* error of closest points grows like
`ε_machine/ε`, while the *distance* between the lines stays accurate. These
behaviours are asserted in `geometry_tests.cpp`.

## 3. Floating-point practices used

* `std::hypot` for norms: no premature overflow/underflow for large/small
  components.
* `atan2(|a×b|, a·b)` for angles: accurate near 0 and π.
* Normalisation checks the length against `zeroLength` *before* dividing.
* All public predicates reject non-finite input instead of letting `NaN`
  comparisons silently evaluate to `false`/`true`.
* No explicit inverses; the only linear solve is an LDLT factorisation of an SPD
  matrix, and its failure is a handled status.
* Finite-difference steps use the representable difference `(x+h) − x`.

## 4. Local degrees of freedom from the Jacobian rank

At a configuration `x` with Jacobian `J` (m × n),

```
DOF_local ≈ n − rank(J) ,        redundant rows = m − rank(J) .
```

*Why this works:* to first order a displacement `δ` preserves all constraints
iff `Jδ = 0`. The solution set of that linear system is the null space of `J`,
whose dimension is `n − rank(J)` (rank–nullity). Those are the directions
along which the points can move without violating the constraints to first
order. `m − rank(J)` counts residual rows that are linear combinations of the
others at `x` (e.g. a duplicated constraint, or the third component of
`û×v̂` at a satisfied configuration).

The rank is the number of singular values `σᵢ > rankTolerance·σ_max` from
`Eigen::BDCSVD`. Default `rankTolerance = 1e-8`:

* the finite-difference Jacobian carries noise of ≈1e-10 relative to its
  entries, so the threshold must sit well above that (otherwise noise creates
  spurious rank);
* it must sit well below any genuine singular value. A *real* but very small
  singular value (< 1e-8 σ_max) — an extremely ill-conditioned yet independent
  constraint — will be reported as rank-deficient. This is inherent to numerical
  rank and is why the threshold is exposed as an option.

### What the estimate does NOT establish

This is a **local, first-order, numerical estimate**. In particular it does not:

* prove global uniqueness: the bracket has DOF = 0 yet two solutions;
* see second-order constraints: with `C = (50,0,0)` at distance 50 from both
  anchors, the Jacobian drops rank (DOF = 1 reported) although `C` is
  geometrically unique, because the free direction only moves the distances at
  second order (`√(50² + y²) − 50 ≈ y²/100`). Demonstrated in
  `solver_tests.cpp` ("tangent configuration");
* classify under/over-constrained systems globally. Whether the constraints
  are *consistent* is answered by the solver status (residual at the end), not by
  the rank; consistent redundancy (duplicate constraints) converges, inconsistent
  constraints end in `Stagnated`;
* apply away from the configuration it was computed at (it is not a property of
  the constraint graph).

Use it as a diagnostic ("this solution still has k free local motions"), not as
a proof of well-constrainedness.

## 5. Tolerance policies and what they cost

| Quantity | Policy | Limitation |
|---|---|---|
| geometry comparisons | absolute + relative (`Tolerances`) | loose in absolute terms at huge magnitudes |
| solver success | **absolute** `max|r| ≤ 1e-9` | at very large scales the floor of double precision (≈ε·|coordinates|) can exceed it, in which case the solver stops as `Stagnated`/`MaxIterations` instead of claiming success; use a larger `residualTolerance`. At the tested scales (1e-3 … 1e6 × the bracket) the default was reachable. |
| step size | absolute + relative | — |
| gradient | absolute | sensitive to residual and Jacobian scaling |
| rank | relative to `σ_max` | see §4 |

`numerical_tests.cpp` asserts the invariant that matters most:
`Converged ⇒ max|r| ≤ residualTolerance`, including for tolerances far below
double precision.

## 6. Where the solver is numerically weak

* **Non-differentiable residual:** `|p−q| − d` at `p = q` has no gradient;
  central differences return zero there (see `limitations.md`).
* **Rank-deficient Jacobians near the solution** (tangent/collinear
  configurations, angle constraints at 0 or π) turn quadratic convergence into
  linear convergence or stagnation; damping keeps the step finite but cannot
  restore speed.
* **Exact symmetry:** if the start lies on a symmetry plane of the problem the
  gradient has no component off the plane and the iteration cannot leave it
  (the bracket started with `y = 0` ends `Stagnated` at a non-solution).
* **Squared condition number** from the normal equations.
