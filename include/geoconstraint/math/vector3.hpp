#pragma once

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>

#include "geoconstraint/numerical/errors.hpp"
#include "geoconstraint/numerical/tolerances.hpp"

namespace geoconstraint {

/// A free vector (direction + magnitude) in R^3. Points are a separate type
/// (`Point3`); the difference of two points is a `Vector3`.
struct Vector3 {
  double x{0.0};
  double y{0.0};
  double z{0.0};

  constexpr Vector3() = default;
  constexpr Vector3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  constexpr Vector3& operator+=(const Vector3& o) {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
  constexpr Vector3& operator-=(const Vector3& o) {
    x -= o.x;
    y -= o.y;
    z -= o.z;
    return *this;
  }
  constexpr Vector3& operator*=(double s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }
  /// Throws InvalidArgumentError when `s` is zero or not finite.
  Vector3& operator/=(double s) {
    if (s == 0.0 || !std::isfinite(s)) {
      throw InvalidArgumentError("Vector3: division by zero or non-finite scalar");
    }
    x /= s;
    y /= s;
    z /= s;
    return *this;
  }

  [[nodiscard]] bool isFinite() const noexcept {
    return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
  }
};

[[nodiscard]] constexpr Vector3 operator+(Vector3 a, const Vector3& b) { return a += b; }
[[nodiscard]] constexpr Vector3 operator-(Vector3 a, const Vector3& b) { return a -= b; }
[[nodiscard]] constexpr Vector3 operator-(const Vector3& a) { return {-a.x, -a.y, -a.z}; }
[[nodiscard]] constexpr Vector3 operator*(Vector3 a, double s) { return a *= s; }
[[nodiscard]] constexpr Vector3 operator*(double s, Vector3 a) { return a *= s; }
[[nodiscard]] inline Vector3 operator/(Vector3 a, double s) { return a /= s; }

[[nodiscard]] constexpr bool operator==(const Vector3& a, const Vector3& b) {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}

[[nodiscard]] constexpr double dot(const Vector3& a, const Vector3& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

/// Right-handed cross product: |a x b| = |a||b| sin(theta), orthogonal to both.
[[nodiscard]] constexpr Vector3 cross(const Vector3& a, const Vector3& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

[[nodiscard]] constexpr double squaredNorm(const Vector3& a) { return dot(a, a); }

/// Euclidean magnitude (computed with hypot, so it neither overflows nor
/// underflows prematurely).
[[nodiscard]] inline double norm(const Vector3& a) { return std::hypot(a.x, a.y, a.z); }

/// Unit vector in the direction of `a`, or std::nullopt if `a` is non-finite or
/// not longer than `tol.zeroLength`.
[[nodiscard]] inline std::optional<Vector3> tryNormalized(
    const Vector3& a, const Tolerances& tol = kDefaultTolerances) {
  if (!a.isFinite()) return std::nullopt;
  const double n = norm(a);
  if (!(n > tol.zeroLength)) return std::nullopt;
  return Vector3{a.x / n, a.y / n, a.z / n};
}

/// Unit vector in the direction of `a`. Throws DegenerateGeometryError for a
/// zero (or non-finite) vector - it never returns an arbitrary direction.
[[nodiscard]] inline Vector3 normalized(const Vector3& a,
                                        const Tolerances& tol = kDefaultTolerances) {
  if (auto n = tryNormalized(a, tol)) return *n;
  throw DegenerateGeometryError("normalized: zero-length or non-finite vector");
}

namespace detail {
/// Throws unless both vectors are finite and longer than tol.zeroLength.
inline void requireDirections(const Vector3& a, const Vector3& b, const Tolerances& tol,
                              const char* what) {
  if (!a.isFinite() || !b.isFinite() || !(norm(a) > tol.zeroLength) ||
      !(norm(b) > tol.zeroLength)) {
    throw DegenerateGeometryError(std::string(what) + ": zero-length or non-finite vector");
  }
}
}  // namespace detail

/// Unsigned angle in [0, pi] between two vectors, computed as
/// atan2(|a x b|, a . b), which stays accurate near 0 and pi (unlike acos).
/// Throws DegenerateGeometryError if either vector is zero.
[[nodiscard]] inline double angleBetween(const Vector3& a, const Vector3& b,
                                         const Tolerances& tol = kDefaultTolerances) {
  detail::requireDirections(a, b, tol, "angleBetween");
  return std::atan2(norm(cross(a, b)), dot(a, b));
}

/// True if |a x b| <= tol.angular * |a||b| (parallel *or* antiparallel).
/// Throws DegenerateGeometryError if either vector is zero: a zero vector is
/// never reported as parallel to anything.
[[nodiscard]] inline bool isParallel(const Vector3& a, const Vector3& b,
                                     const Tolerances& tol = kDefaultTolerances) {
  detail::requireDirections(a, b, tol, "isParallel");
  return norm(cross(a, b)) <= tol.angular * norm(a) * norm(b);
}

/// True if |a . b| <= tol.angular * |a||b|. Throws on zero vectors.
[[nodiscard]] inline bool isPerpendicular(const Vector3& a, const Vector3& b,
                                          const Tolerances& tol = kDefaultTolerances) {
  detail::requireDirections(a, b, tol, "isPerpendicular");
  return std::abs(dot(a, b)) <= tol.angular * norm(a) * norm(b);
}

/// |a-b| <= length + relative * max(|a|,|b|).
[[nodiscard]] inline bool approxEqual(const Vector3& a, const Vector3& b,
                                      const Tolerances& tol = kDefaultTolerances) {
  if (!a.isFinite() || !b.isFinite()) return false;
  return norm(a - b) <= lengthTolerance(tol, std::max(norm(a), norm(b)));
}

}  // namespace geoconstraint
