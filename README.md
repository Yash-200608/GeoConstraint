# GeoConstraint

A small, tested C++20 library for **3D computational geometry** and **geometric
constraint solving**: vectors, points, rotations and rigid transforms; lines,
planes, projections, distances and intersections; seven geometric constraint
types; and a Levenberg–Marquardt least-squares solver with structured
convergence, residual and degree-of-freedom diagnostics.

> **Scope.** GeoConstraint is an independent, educational engine. It is not a
> production CAD constraint solver, contains no third-party proprietary code, and
> makes no claim of compatibility with any commercial component. See
> [`docs/limitations.md`](docs/limitations.md).

## Why this problem matters in engineering geometry

Parametric CAD describes a part by *relationships* — "this hole is 12 mm from
that edge", "these faces are perpendicular" — rather than by raw coordinates.
Changing one dimension means finding new positions that satisfy *all*
relationships at once: a system of nonlinear equations that may be
under-constrained (parts can still move), over-constrained or inconsistent, have
several valid solutions, or be numerically delicate near degenerate geometry.
A geometric constraint solver is the component that does this. This project
implements the essential pieces end to end — geometry kernel, constraint
residuals, a damped Gauss–Newton solver, and honest diagnostics — at a scale
small enough to read and verify.

## What is supported

### Geometry (all on **infinite** lines and planes)

| Operation | API |
|---|---|
| Vector algebra: `+ − * /`, dot, cross, norm, safe normalisation, unsigned angle, parallel/perpendicular tests, approximate equality | `Vector3` |
| Points (distinct from vectors), distance | `Point3` |
| 3×3 matrices, axis–angle rotations (Rodrigues) | `Matrix3` |
| Rigid transforms: apply to points and vectors, compose, invert, rotate about a point | `RigidTransform3` |
| Lines `Line3`, planes `Plane3` (validated factories, unit direction/normal) | `geometry/line3.hpp`, `plane3.hpp` |
| Projection of a point onto a line / plane | `projectPointOnLine`, `projectPointOnPlane` |
| Point–line, point–plane (unsigned/signed), point–point distance | `distances.hpp` |
| Line–plane intersection: unique point / none / line in plane / invalid | `intersectLinePlane` |
| Closest points and distance of two lines: skew / intersecting / parallel / coincident / invalid | `closestPointsBetweenLines` |

Rays and segments are **not** implemented; nothing silently clamps to a segment.
Zero vectors, non-finite values, collinear plane definitions, etc. throw
`DegenerateGeometryError` / `InvalidArgumentError`; they are never "repaired".

### Constraints

Points are referenced by stable ids (`PointId`). Residual vectors are
concatenated in insertion order.

| Constraint | Residual | # |
|---|---|---|
| `CoincidenceConstraint(p,q)` | `p − q` | 3 |
| `FixedPointConstraint(p,p₀)` | `p − p₀` | 3 |
| `DistanceConstraint(p,q,d)` | `|p − q| − d` | 1 |
| `PointOnPlaneConstraint(p,plane)` | `n·(p − p₀)` | 1 |
| `ParallelConstraint(u,v)` | `û × v̂` (parallel *and* antiparallel satisfy it) | 3 |
| `PerpendicularConstraint(u,v)` | `û · v̂` | 1 |
| `AngleConstraint(u,v,θ)` | `û · v̂ − cos θ`, unsigned `θ ∈ [0, π]` of *directed* directions | 1 |

Directions are `DirectionRef{from,to}` (vector `P[to] − P[from]`). A zero-length
direction is an error, not a trivially satisfied constraint. Full definitions and
conventions: [`docs/mathematical_foundations.md`](docs/mathematical_foundations.md).

### Solver

Levenberg–Marquardt on `½‖r(x)‖²` with a finite-difference Jacobian
(forward or central), the damped normal equations `(JᵀJ + λI)Δx = −Jᵀr`
solved with `Eigen::LDLT` (no explicit inverse), gain-ratio step acceptance, and
adaptive damping. Outcomes: `Converged` (only if `max|rᵢ| ≤ residualTolerance`),
`MaxIterations`, `Stagnated`, `InvalidInput`, `NumericalFailure`. Every result
carries per-constraint residuals, counts, a cost history and an SVD-based
**local** estimate of the Jacobian rank and degrees of freedom.
Details: [`docs/solver_design.md`](docs/solver_design.md).

## Usage

```cpp
#include "geoconstraint/constraints/distance_constraint.hpp"
#include "geoconstraint/constraints/fixed_point_constraint.hpp"
#include "geoconstraint/constraints/point_on_plane_constraint.hpp"
#include "geoconstraint/solver/solver.hpp"
using namespace geoconstraint;

Problem problem;
PointId A = problem.addPoint({0, 0, 0}, "A");
PointId B = problem.addPoint({100, 0, 0}, "B");
PointId C = problem.addPoint({40, 30, 10}, "C");       // initial guess

problem.add<FixedPointConstraint>(A, Point3{0, 0, 0});  // anchors
problem.add<FixedPointConstraint>(B, Point3{100, 0, 0});
problem.add<DistanceConstraint>(C, A, std::sqrt(5000.0));
problem.add<DistanceConstraint>(C, B, std::sqrt(5000.0));
problem.add<PointOnPlaneConstraint>(C, Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1}));

SolverResult result = solve(problem);                   // SolverOptions{} defaults
if (result.converged()) {
  Point3 c = result.points[C.index];                    // (50, 50, 0)
}
std::cout << formatReport(problem, result);
```

CMake consumers: `add_subdirectory(GeoConstraint)` then
`target_link_libraries(app PRIVATE GeoConstraint::geoconstraint)`.

## Build and test

Requirements: a C++20 compiler with `<span>`, `<compare>` and `<numbers>` (expected:
GCC ≥ 11, Clang ≥ 14, MSVC 2019 16.11+/2022; only Clang 21 was actually tried),
CMake ≥ 3.21. **Eigen 3.4** and **Catch2 3.7** are used from a system install if
`find_package` locates them; otherwise CMake downloads the pinned release
archives (SHA-256 verified), so the first configure needs network access.

Windows PowerShell (Visual Studio generator, multi-config):

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Linux (single-config generator):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

CMake options: `GEOCONSTRAINT_BUILD_TESTS` (ON), `GEOCONSTRAINT_BUILD_EXAMPLES`
(ON), `GEOCONSTRAINT_WARNINGS_AS_ERRORS` (OFF). With a single-config generator
`CMAKE_BUILD_TYPE` is honoured; with Visual Studio it is ignored in favour of
`--config`.

### What has actually been verified

The project was built and tested on Windows 11 with CMake 4.4.4, Ninja and the
Clang 21 front end shipped with Zig 0.16 (`zig c++`), configured with
`-DCMAKE_BUILD_TYPE=Release -DGEOCONSTRAINT_WARNINGS_AS_ERRORS=ON`:

* Clean build: no warnings with `-Wall -Wextra -Wpedantic -Wshadow -Wconversion
  -Wsign-conversion -Wnon-virtual-dtor -Woverloaded-virtual -Wdouble-promotion`.
* `ctest`: **81/81 passed** (78 Catch2 test cases + 3 example programs).
* A Debug build with `-fsanitize=undefined` also passed all 81 with no reports.
* Three deliberate defects injected into the implementation (distance residual
  sign, line–line parameter sign, LM predicted-reduction sign) were each caught
  by the suite (14, 3 and 16 failing tests), then reverted.

**Not verified locally:** MSVC and GCC builds, Linux. They are covered by
`.github/workflows/ci.yml` (Linux GCC, Linux Clang with `-Werror`, Windows MSVC),
which has not been run yet. The first CI run on those compilers may surface
warnings or portability issues the author could not see.

## Examples

All examples are built by default, registered in CTest, and exit non-zero if
their own checks fail.

| Program | Shows |
|---|---|
| `geometry_demo` | projections, distances, intersections, closest lines, transforms, each checked against the known answer |
| `bracket_assembly [x y z]` | the end-to-end demonstration below; optional custom starting guess for C |
| `solver_diagnostics` | underdetermined (DOF = 2), inconsistent, iteration-limited and invalid inputs |

### Bracket assembly

`A = (0,0,0)` and `B = (100,0,0)` are pinned by fixed-point constraints. `C` is
free, constrained by `|CA| = |CB| = √5000` and `z_C = 0`, which has two
solutions, `C = (50, ±50, 0)`. Which one is found depends only on the start
(`y > 0` → `+50`, `y < 0` → `−50`). Output of `bracket_assembly` (first case,
abridged; counts and last digits can vary slightly between compilers):

```
=== upper branch (guess y > 0) ===
Initial guess for C : (40.000000, 30.000000, 10.000000)
Status        : Converged - all residuals within residualTolerance
Iterations    : 5 accepted, 0 rejected steps, 114 residual evaluations
Residual      : |r|_2 = 2.329e-13 (initial 2.230e+01), |r|_inf = 1.421e-13
Points:
  A = (0.000000000, 0.000000000, -0.000000000)
  B = (100.000000000, 0.000000000, -0.000000000)
  C = (50.000000000, 50.000000000, 0.000000000)
Jacobian (local numerical estimate): 9 x 9, rank = 9, local DOF = 0, ...
Anchor drift: A 4.547e-14, B 1.961e-13
CHECK PASSED
```

Run it: `build/examples/bracket_assembly` (Linux/Ninja) or
`build\examples\Release\bracket_assembly.exe` (Visual Studio).

## Numerical conventions

* All thresholds are in one struct, `Tolerances`: `zeroLength = 1e-12`,
  `length = 1e-9`, `relative = 1e-9`, `angular = 1e-9`. Length comparisons use
  `length + relative · scale`. Intended coordinate scale ≈ 1e-3 … 1e6.
* Solver success is **absolute**: `max|rᵢ| ≤ 1e-9` by default
  (`SolverOptions::residualTolerance`).
* Angles are unsigned in `[0, π]`, computed with `atan2(|a×b|, a·b)`.
* The local DOF estimate `n − rank(J)` uses `rank = #{σᵢ > 1e-8 σ_max}` and is a
  first-order, local estimate only.

See [`docs/numerical_robustness.md`](docs/numerical_robustness.md).

## Repository layout

```
include/geoconstraint/
  numerical/   tolerances.hpp, errors.hpp
  math/        vector3, point3, matrix3, rigid_transform3          (header-only)
  geometry/    line3, plane3, projections, distances, intersections
  constraints/ constraint.hpp + one header per constraint
  solver/      problem, solver, solver_options, solver_result, diagnostics
src/           geometry/  constraints/  solver/        (static library)
tests/         math, geometry, constraint, solver, numerical (Catch2)
examples/      geometry_demo, bracket_assembly, solver_diagnostics
docs/          mathematical_foundations, architecture, solver_design,
               numerical_robustness, limitations
.github/workflows/ci.yml
```

`docs/architecture.md` explains the layering and error-handling policy. (The
suggested `src/cli/` was omitted: `bracket_assembly [x y z]` serves that role.)

## Known limitations (summary)

* Local solver: finds *a* solution near the start; which branch is found
  depends on the start. Exactly symmetric starts (e.g. the bracket with `y = 0`)
  are stationary points and end as `Stagnated`.
* Only points are unknowns; planes are constants. Seven constraint types.
* Dense finite-difference Jacobian: suited to small systems.
* `Stagnated` does not prove infeasibility; the DOF estimate is local and
  first-order (it reports 1 free DOF for the tangent bracket configuration
  although the point is geometrically unique).
* Distance residuals are non-differentiable at coincident points; angle
  constraints at θ = 0 or π are first-order singular.
* No sparsity, no analytic Jacobians, no conflict identification.

Full list: [`docs/limitations.md`](docs/limitations.md).

## Possible future extensions

Analytic Jacobians and sparse/decomposed solving; QR-based LM step and
Marquardt scaling; more entities (lines, circles, planes as unknowns) and
constraints (tangency, point-on-line, symmetry); signed angles; structural
DOF analysis and conflict detection; rays/segments; installable CMake package
and Conan/vcpkg recipes.

## License

MIT — see [`LICENSE`](LICENSE).
