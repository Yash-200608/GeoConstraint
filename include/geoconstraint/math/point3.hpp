#pragma once

#include <algorithm>
#include <cmath>

#include "geoconstraint/math/vector3.hpp"

namespace geoconstraint {

/// A location in R^3. Distinct from `Vector3`: points cannot be added or
/// scaled; only `Point - Point -> Vector` and `Point +/- Vector -> Point` exist.
struct Point3 {
  double x{0.0};
  double y{0.0};
  double z{0.0};

  constexpr Point3() = default;
  constexpr Point3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  constexpr Point3& operator+=(const Vector3& v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }
  constexpr Point3& operator-=(const Vector3& v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  }

  /// Position vector from the origin to this point.
  [[nodiscard]] constexpr Vector3 positionVector() const { return {x, y, z}; }

  [[nodiscard]] bool isFinite() const noexcept {
    return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
  }
};

[[nodiscard]] constexpr Point3 operator+(Point3 p, const Vector3& v) { return p += v; }
[[nodiscard]] constexpr Point3 operator-(Point3 p, const Vector3& v) { return p -= v; }
[[nodiscard]] constexpr Vector3 operator-(const Point3& a, const Point3& b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
[[nodiscard]] constexpr bool operator==(const Point3& a, const Point3& b) {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}

[[nodiscard]] constexpr double squaredDistance(const Point3& a, const Point3& b) {
  return squaredNorm(a - b);
}

/// Euclidean distance between two points.
[[nodiscard]] inline double distance(const Point3& a, const Point3& b) { return norm(a - b); }

/// Points are equal if their distance is within the effective length tolerance.
[[nodiscard]] inline bool approxEqual(const Point3& a, const Point3& b,
                                      const Tolerances& tol = kDefaultTolerances) {
  if (!a.isFinite() || !b.isFinite()) return false;
  return distance(a, b) <=
         lengthTolerance(tol, std::max(norm(a.positionVector()), norm(b.positionVector())));
}

}  // namespace geoconstraint
