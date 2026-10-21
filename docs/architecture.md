# Architecture

GeoConstraint is a small layered C++20 library. Dependencies point downward only.

```
 examples/ tests/                      (consumers)
       │
 solver/        Problem, solve(), SolverOptions, SolverResult, diagnostics
       │
 constraints/   Constraint interface + 7 concrete constraints
       │
 geometry/      Line3, Plane3, projections, distances, intersections
       │
 math/          Vector3, Point3, Matrix3, RigidTransform3
       │
 numerical/     Tolerances, error types
```

## Layers

### `numerical/`
* `tolerances.hpp` — the single `Tolerances` policy (`constexpr` default), plus
  `approxEqual(double,…)` and `lengthTolerance`.
* `errors.hpp` — `GeoConstraintError` base; `InvalidArgumentError` (value out of
  domain) and `DegenerateGeometryError` (valid numbers, unusable geometry).

### `math/` (header-only)
Value types with no allocation. `Point3` and `Vector3` are deliberately
unconvertible; `Matrix3` is a row-major 3×3 with Rodrigues rotations and
`isRotation`; `RigidTransform3` validates its rotation and stores `R, t`.

### `geometry/`
`Line3` and `Plane3` are immutable classes created only through validating
factories — an existing object is never degenerate. Free functions implement
projections (header-only), distances and line–line closest points
(`distances.cpp`) and line–plane intersection (`line_plane.cpp`). Results with
several possible outcomes are returned as an enum plus `std::optional` payloads
(`LinePlaneIntersection`, `LineLineClosest`) so the "no unique answer" cases
cannot be mistaken for a point.

### `constraints/`
`Constraint` is an abstract interface:

```cpp
typeName();  residualCount();  referencedPoints();
evaluate(std::span<const Point3> points, std::span<double> residuals);
```

Constraints hold **point ids** (`PointId`, an index into the problem's point
list) and constants, never coordinates; the solver owns the coordinates. This is
why a coordinate has exactly one definition and why moving a point updates every
constraint that references it. Direction-based constraints take
`DirectionRef{from,to}` (two ids). `evaluate` is const and side-effect free, so
constraints are evaluated independently and their residuals concatenated in a
deterministic order. Concrete classes are in one header each; their
implementations share `src/constraints/constraints.cpp`.

### `solver/`
* `Problem` — owns initial positions, names, and `unique_ptr<Constraint>`s
  (RAII; move-only). Adding a constraint validates its point ids. Provides the
  state codec (`encode`/`decode`) and the combined residual `evaluateResiduals`.
* `solve` (`solver.cpp`) — Levenberg–Marquardt; never throws for bad problems,
  returns a `SolverResult`.
* `diagnostics` — `ConstraintResidual`, `SystemAnalysis`, SVD rank analysis,
  and the human-readable `formatReport`.
* `finiteDifferenceJacobian` and `analyzeSystem` are public for testing and for
  analysing arbitrary configurations.

## Error-handling policy

| Layer | Mechanism |
|---|---|
| math / geometry / constraint *construction* | throws (`InvalidArgumentError`, `DegenerateGeometryError`) — a programming or data error the caller must handle |
| geometry queries with several legitimate outcomes | enum classification (`LinePlaneRelation`, `LineLineRelation`) |
| solver | **no exceptions for bad problems**: structured `SolverStatus` + `message`; exceptions from constraints are caught at the solver boundary and become `InvalidInput` (initial state) or rejected trial steps / `NumericalFailure` (mid-solve) |

Programming errors in the API (null constraint, out-of-range id) throw from
`Problem::addConstraint` immediately.

## Ownership and lifetimes

* `Problem` owns its constraints; `solve` takes `const Problem&` and returns an
  independent `SolverResult` (it stores no references to the problem).
* `PointId` is a plain index: valid for the `Problem` that issued it. Ids are
  stable because points are never removed.
* No global mutable state; `kDefaultTolerances` is `constexpr`.

## Build structure

One static library `geoconstraint` (sources in `src/`), header-only math layer,
public include directory `include/`, Eigen as a public dependency (it appears in
solver headers). Tests build into one Catch2 executable registered with CTest via
`catch_discover_tests`; each example is registered as a CTest test and exits
non-zero if its own self-checks fail. Deviations from the suggested layout: the
`src/cli/` directory is omitted (the `bracket_assembly` example accepts a custom
starting guess and takes its role), and small related implementations share a
`.cpp` (all constraints in `constraints.cpp`, line–plane and distances in two
geometry files).
