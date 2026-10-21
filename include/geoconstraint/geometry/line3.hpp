#pragma once

#include "geoconstraint/math/point3.hpp"

namespace geoconstraint {

/// An INFINITE line {origin + t * direction : t in R} with a unit direction.
/// (Not a ray, not a segment: every operation in this library that takes a
/// `Line3` treats it as unbounded in both directions.)
///
/// Invariant: `direction()` is finite with length 1 (to rounding); `origin()`
/// is finite. Only the named factories can create a line, so a `Line3` is
/// never degenerate.
class Line3 {
 public:
  /// Throws DegenerateGeometryError if `direction` is zero/non-finite and
  /// InvalidArgumentError if `origin` is non-finite.
  [[nodiscard]] static Line3 fromPointDirection(const Point3& origin, const Vector3& direction,
                                                const Tolerances& tol = kDefaultTolerances) {
    if (!origin.isFinite()) throw InvalidArgumentError("Line3: non-finite origin");
    return Line3(origin, normalized(direction, tol));
  }

  /// Line through two distinct points (direction b - a). Throws
  /// DegenerateGeometryError if the points coincide (within tol.zeroLength).
  [[nodiscard]] static Line3 fromTwoPoints(const Point3& a, const Point3& b,
                                           const Tolerances& tol = kDefaultTolerances) {
    return fromPointDirection(a, b - a, tol);
  }

  [[nodiscard]] const Point3& origin() const noexcept { return origin_; }
  [[nodiscard]] const Vector3& direction() const noexcept { return direction_; }

  /// Point at signed arc-length parameter t (the direction is a unit vector).
  [[nodiscard]] Point3 pointAt(double t) const { return origin_ + direction_ * t; }

  /// Signed parameter of the orthogonal projection of `p` onto the line.
  [[nodiscard]] double parameterOf(const Point3& p) const { return dot(p - origin_, direction_); }

 private:
  Line3(const Point3& o, const Vector3& d) : origin_(o), direction_(d) {}
  Point3 origin_;
  Vector3 direction_;
};

}  // namespace geoconstraint
