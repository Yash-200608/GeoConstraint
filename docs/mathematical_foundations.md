# Mathematical foundations

This document derives every formula the library uses and explains *why* it is
correct. Notation: bold-free plain symbols; `a·b` is the dot product, `a×b` the
cross product, `|a|` the Euclidean norm, and a hat (`û`) a unit vector.

Points (`Point3`) and vectors (`Vector3`) are different types. Only
`point − point → vector` and `point ± vector → point` exist, mirroring the
affine-geometry rule that positions cannot be added or scaled while
displacements can.

## 1. Dot and cross products

### Dot product

`a·b = ax·bx + ay·by + az·bz`. Geometrically `a·b = |a||b|cos θ`, where `θ` is
the angle between the vectors.

*Why:* apply the law of cosines to the triangle with sides `a`, `b`, `a−b`:
`|a−b|² = |a|² + |b|² − 2|a||b|cos θ`. Expanding `|a−b|² = |a|² + |b|² − 2 a·b`
and cancelling gives the identity. Consequences used by the library:

* `a·b = 0` ⇔ perpendicular (for nonzero vectors). `isPerpendicular` tests
  `|a·b| ≤ tol·|a||b|`, i.e. `|cos θ| ≤ tol`, which is scale invariant.
* For a unit vector `d`, `(p−o)·d` is the signed length of the component of
  `p−o` along `d` (used by projections).

### Cross product

`a×b = (ay bz − az by, az bx − ax bz, ax by − ay bx)`. It is orthogonal to both
inputs, right-handed, and `|a×b| = |a||b| sin θ`.

*Why orthogonal:* `a·(a×b)` expands to six terms that cancel in pairs.
*Why the magnitude:* the Lagrange identity
`|a×b|² + (a·b)² = |a|²|b|²` (checked by `math_tests.cpp`) combined with
`(a·b)² = |a|²|b|² cos²θ` gives `|a×b|² = |a|²|b|² sin²θ`.

Consequences:

* `a×b = 0` ⇔ parallel *or antiparallel* (for nonzero vectors).
  `isParallel` tests `|a×b| ≤ tol·|a||b|`, i.e. `sin θ ≤ tol`.
* The unsigned angle is computed as `θ = atan2(|a×b|, a·b)`. Unlike
  `acos(â·b̂)`, whose argument is flat near ±1 (`acos(1−δ) ≈ √(2δ)` amplifies
  rounding of order 1e-16 into 1e-8 errors for small angles), `atan2` uses both
  the sine and cosine and is accurate for the whole range `[0, π]`. This is
  verified by `numerical_tests.cpp` at θ = 1e-8.

### Rotation matrix (Rodrigues)

For a unit axis `u` and angle `a`:
`R = I + sin a·K + (1 − cos a)·K²`, `K = [u]ₓ` (so `K v = u×v`). A matrix is a
proper rotation if `RᵀR = I` and `det R = +1`; `RigidTransform3` refuses
anything else (a reflection has `det = −1`, a scaling breaks `RᵀR = I`).

### Rigid transforms

`x ↦ R x + t`. Composition `(A∘B)(x) = A(B(x)) = R_A R_B x + (R_A t_B + t_A)`;
the library writes `A * B` for "apply `B` first". The inverse is
`x ↦ Rᵀ(x − t)`. Directions transform with `R` only, because the translation
cancels in a difference of points: `(Rp+t) − (Rq+t) = R(p−q)`.

## 2. Projections, distances and intersections

All objects are **infinite**: `Line3` is `{o + t d : t ∈ ℝ}` (not a ray, not a
segment) and `Plane3` is `{p : n·(p − o) = 0}`; `d` and `n` are unit vectors.
Nothing clamps to a segment; a projection can lie arbitrarily far from the
defining points.

### Point onto a line

Minimise `f(t) = |p − (o + t d)|²`. `f'(t) = −2 d·(p − o − t d) = 0` and
`d·d = 1` give `t* = d·(p − o)`, so the foot point is
`o + ((p−o)·d) d`. The residual `p − foot` is orthogonal to `d` by construction.
Distance: `|p − foot| = |(p−o) − ((p−o)·d)d| = |(p−o)×d|` (the component of
`p−o` orthogonal to `d`; its norm is `|p−o| sin θ`).

### Point onto a plane

The signed distance is `s = n·(p − o)`. Moving `p` by `−s n` removes exactly
the normal component, so the projection is `p − s n` and the distance is `|s|`.
Check: `n·(p − s n − o) = s − s(n·n) = 0`.

### Line–plane intersection

Insert the line into the plane equation: `n·(o_ℓ + t d − o_π) = 0`, so
`t = n·(o_π − o_ℓ) / (n·d)`.

Because both `n` and `d` are unit vectors, `n·d` is the sine of the angle
between the line and the plane. The classification is therefore:

| `|n·d|` | Offset `n·(o_ℓ − o_π)` | Result |
|---|---|---|
| `> tol.angular` | — | `UniquePoint` at the parameter above |
| `≤ tol.angular` | within length tolerance | `LineInPlane` (every point is common) |
| `≤ tol.angular` | otherwise | `NoIntersection` (parallel, off the plane) |
| non-finite anywhere | — | `Invalid` |

`point`/`parameter` are `std::optional` and are only present for a unique
intersection: no representative point is invented for the non-unique cases.

### Closest points of two lines

Lines `A: oA + s dA`, `B: oB + t dB`, `w = oB − oA`, `n = dA×dB`.
The segment joining the closest points is the common perpendicular, hence
parallel to `n`:

```
s dA − t dB = w + λ n .
```

Cross both sides with `dB` and dot with `n`; since `(n×dB)·n = 0` and
`dB×dB = 0`:

```
s (dA×dB)·n = (w×dB)·n   ⇒   s = ((w×dB)·n) / |n|² .
```

Crossing with `dA` instead gives `t = ((w×dA)·n) / |n|²`. The distance is the
projection of `w` on the unit common normal: `|w·n| / |n|`.

* `|n| = sin(angle)`. If `|n| ≤ tol.angular` the lines are parallel and the
  formulas above would divide by ~0, so the parallel branch is used instead:
  the distance is that of `oA` to line B, and `Coincident` is reported when it
  is within the length tolerance, `Parallel` otherwise. The closest pair is not
  unique there; a representative pair is returned and documented as such.
* Otherwise `Intersecting` if the distance is within tolerance, else `Skew`.

## 3. Constraint residual definitions

Every constraint is satisfied exactly when all its scalar residuals are zero.
The residual count is a property of the *type*, never of the geometry, so the
global residual vector has a fixed layout.

| Constraint | Residual | Count | Units |
|---|---|---|---|
| `CoincidenceConstraint(p,q)` | `r = p − q` | 3 | length |
| `FixedPointConstraint(p,p₀)` | `r = p − p₀` | 3 | length |
| `DistanceConstraint(p,q,d)` | `r = |p − q| − d` | 1 | length |
| `PointOnPlaneConstraint(p,Π)` | `r = n·(p − o)` | 1 | length (signed distance) |
| `ParallelConstraint(u,v)` | `r = û × v̂` | 3 | dimensionless (`sin θ`) |
| `PerpendicularConstraint(u,v)` | `r = û · v̂` | 1 | dimensionless (`cos θ`) |
| `AngleConstraint(u,v,θ)` | `r = û · v̂ − cos θ` | 1 | dimensionless |

Directions `u`, `v` are `DirectionRef{from,to}` pairs of point ids, with vector
`P[to] − P[from]` normalised before use. Using unit directions makes direction
residuals independent of segment length, so they are comparable in size with
each other and bounded by 1.

### Direction conventions

* **Directed:** reversing the order of the two points of a `DirectionRef`
  flips its direction. This maps an angle `θ` to `π − θ`.
* **Angle range and sign:** `θ ∈ [0, π]`, unsigned; the constraint is
  symmetric in its two directions. `θ = 0` means same direction, `θ = π`
  antiparallel, `θ = π/2` perpendicular.
* **Antiparallel directions:** satisfy `ParallelConstraint` (`û×v̂ = 0`) and
  `AngleConstraint(π)`, but *not* `AngleConstraint(0)`, whose residual is `−2`.
* **Zero directions are errors.** A direction shorter than
  `Tolerances::zeroLength` (1e-12) or containing a non-finite coordinate raises
  `DegenerateGeometryError`. Silently normalising it would produce `NaN`, and a
  zero vector would trivially satisfy `u×v = 0` and `u·v = 0`.
* **Why cos-based angles are singular at 0 and π:** `∂(cos θ)/∂θ = −sin θ`
  vanishes there, so an angle constraint with `θ ∈ {0, π}` has a rank-deficient
  Jacobian and converges only linearly (see `limitations.md`).
  `ParallelConstraint` has the same first-order structure: at a satisfied
  configuration its 3×n Jacobian has rank 2, not 3 (`û×v̂` varies only in the
  plane orthogonal to `û` there), which the diagnostics report as one redundant
  residual row.

## 4. Optimisation problem

Let `x ∈ ℝⁿ` hold all point coordinates (`n = 3 × #points`, point `k` occupying
entries `3k..3k+2`) and `r(x) ∈ ℝᵐ` the concatenated residuals. The solver
minimises

```
F(x) = ½ |r(x)|² .
```

`F ≥ 0` and `F(x) = 0` ⇔ every constraint is satisfied, so zero-residual
minima are exactly the solutions. Fixed points are *not* removed from `x`:
`FixedPointConstraint` pins them softly, which keeps the formulation uniform
and lets the diagnostics treat an anchor like any other constraint.

The Levenberg–Marquardt update, the finite-difference Jacobian, the stopping
rules and the rank/DOF analysis are derived in [`solver_design.md`](solver_design.md)
and [`numerical_robustness.md`](numerical_robustness.md).
