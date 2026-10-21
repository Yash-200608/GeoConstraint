#pragma once

#include "geoconstraint/math/point3.hpp"

namespace geoconstraint {

/// An INFINITE plane {p : n . (p - origin) = 0} with a unit normal.
///
/// The normal also defines the sign convention of signed distances: points on
/// the side the normal points to have positive signed distance.
class Plane3 {
 public:
  /// Throws DegenerateGeometryError for a zero/non-finite normal and
  /// InvalidArgumentError for a non-finite origin. The normal is normalized.
  [[nodiscard]] static Plane3 fromPointNormal(const Point3& origin, const Vector3& normal,
                                              const Tolerances& tol = kDefaultTolerances) {
    if (!origin.isFinite()) throw InvalidArgumentError("Plane3: non-finite origin");
    return Plane3(origin, normalized(normal, tol));
  }

  /// Plane through three points, normal (b-a) x (c-a) (right-hand rule in the
  /// order a, b, c). Throws DegenerateGeometryError if the points are
  /// coincident or collinear: |(b-a) x (c-a)| <= tol.angular * |b-a| |c-a|.
  [[nodiscard]] static Plane3 fromThreePoints(const Point3& a, const Point3& b, const Point3& c,
                                              const Tolerances& tol = kDefaultTolerances) {
    const Vector3 ab = b - a;
    const Vector3 ac = c - a;
    const Vector3 n = cross(ab, ac);
    if (!n.isFinite() || !(norm(ab) > tol.zeroLength) || !(norm(ac) > tol.zeroLength) ||
        norm(n) <= tol.angular * norm(ab) * norm(ac)) {
      throw DegenerateGeometryError("Plane3::fromThreePoints: collinear or coincident points");
    }
    return fromPointNormal(a, n, tol);
  }

  [[nodiscard]] const Point3& origin() const noexcept { return origin_; }
  [[nodiscard]] const Vector3& normal() const noexcept { return normal_; }

  /// n . (p - origin): positive on the side the normal points to.
  [[nodiscard]] double signedDistance(const Point3& p) const { return dot(normal_, p - origin_); }

 private:
  Plane3(const Point3& o, const Vector3& n) : origin_(o), normal_(n) {}
  Point3 origin_;
  Vector3 normal_;
};

}  // namespace geoconstraint
