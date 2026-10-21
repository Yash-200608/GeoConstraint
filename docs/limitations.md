# Limitations

GeoConstraint is an independent, educational engine. It is **not** a
production CAD constraint solver, makes no claim of compatibility with any
commercial component (including Siemens D-Cubed), and should not be used where
a robust, general solver is required. The items below are known limits of this
release.

## Scope of the model

* **Entities:** points (the only unknowns), and *fixed* planes (as constants of
  `PointOnPlaneConstraint`). There are no free lines, planes, circles, arcs,
  curves or surfaces as solver entities; line directions are derived from two
  points. Lines and planes in the geometry layer are infinite; rays and
  segments are not implemented.
* **Constraints:** exactly the seven listed (coincidence, fixed point, distance,
  point-on-plane, parallel, perpendicular, angle). Missing, among others:
  point-on-line, point-on-circle, tangency, symmetry, equal length, plane–plane
  relations, driven dimensions, weights/priorities.
* **No tolerance/units system:** coordinates are plain doubles in consistent,
  unspecified units; the defaults assume ≈1e-3 … 1e6.
* **Fixed points are soft constraints,** satisfied to the solver tolerance
  (≈1e-9 by default), not removed from the unknowns.

## Constraint-specific caveats

* **Distance:** `|p−q|−d` is not differentiable at `p = q`; central differences
  give a zero gradient there, so a start with coincident points and `d > 0`
  stagnates (or is invalid for related direction constraints).
* **Angle:** angles are unsigned, in `[0, π]`, and cos-based. At `θ = 0` and
  `θ = π` the gradient vanishes, the Jacobian loses rank and convergence from a
  perturbed start is slow/linear. No oriented (signed) angles; no angle between
  a line and a plane.
* **Parallel:** accepts antiparallel directions. 3 residuals, only 2
  independent near a solution — it reports one redundant row by design.
* **Perpendicular, parallel, angle with a degenerate direction** are errors
  (`InvalidInput` when present at the start). If a *trial* iterate makes a
  direction degenerate the step is rejected; if the problem *needs* to pass
  through such a state it can fail.

## Solver

* **Local method:** it finds *a* nearby solution of a nonlinear system, never a
  global one, and it cannot enumerate solutions. Which branch of a multi-solution
  problem is found (bracket: `y > 0` vs `y < 0`) depends on the start.
* **Symmetric stationary points:** a start exactly on a symmetry plane, or any
  other stationary point of `½|r|²` with nonzero residual, cannot be left by a
  gradient-based method. The result is `Stagnated` (verified: bracket with
  `C = (50,0,0)` or `y = 0`). Perturb the start off the symmetry.
* **Stagnation is not infeasibility:** `Stagnated` means "no more progress" —
  inconsistent constraints are the usual cause, but a local minimum of a
  feasible system gives the same status. The library does not prove
  infeasibility or identify *which* constraints conflict (it reports
  per-constraint residuals to help a human find them).
* **Dense finite-difference Jacobian:** `O(n)` residual evaluations per
  iteration and dense `O(n³)` factorisations; fine for dozens of points, not for
  thousands. No sparsity, no analytic Jacobians, no problem decomposition into
  independent sub-systems.
* **Identity damping** is not scale invariant; residuals of very different
  magnitude or custom constraints with huge Jacobian entries may need rescaling.
* **Normal equations** square the condition number of `J`.
* **Absolute residual tolerance:** at very large coordinate scales the default
  `1e-9` may be below the achievable accuracy of double precision (the solver
  then reports non-convergence rather than succeeding falsely); pass a suitable
  `residualTolerance`.
* **Not thread-safe by shared mutation, but re-entrant:** `solve` takes a
  `const Problem&` and keeps no global state, so independent solves can run on
  different threads.

## Diagnostics

* **Local DOF / rank is a first-order numerical estimate** (see
  `numerical_robustness.md` §4): it can overstate freedom at tangent/collinear
  configurations (second-order-constrained motions), cannot detect multiple
  solutions, and depends on `rankTolerance` (genuinely tiny but independent
  singular values are reported as rank loss).
* **No structural/graph analysis:** over-/under-constraint is not classified
  combinatorially; redundancy is only the count `m − rank` at the final point.
* **No identification of conflicting subsets** of constraints.
* **Unsupported:** the DOF estimate is computed for the *final* configuration
  even when the solve fails; in that case it describes the (non-solution) final
  iterate.

## Geometry layer

* Infinite lines/planes only; near-parallel classification is by the angular
  tolerance (1e-9 rad), after which results can contain very distant points.
* Tolerances are global scalars chosen for ≈1e-3 … 1e6 coordinates; there is no
  per-model scaling.
* No polygons, meshes, curves, surfaces, booleans or topology.

## Build and environment

* The test and example suites were developed and run with the Clang-based Zig
  toolchain on Windows (see README for exactly what was verified). MSVC, GCC,
  and Linux builds are configured (CI) but are not verified by the author's
  local run.
