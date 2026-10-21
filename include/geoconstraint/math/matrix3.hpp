#pragma once

#include <array>
#include <cmath>

#include "geoconstraint/math/vector3.hpp"

namespace geoconstraint {

/// Dense 3x3 matrix, row-major. Used for rotations (see `isRotation`).
class Matrix3 {
 public:
  /// Identity matrix.
  constexpr Matrix3() : m_{1, 0, 0, 0, 1, 0, 0, 0, 1} {}
  /// Rows given explicitly.
  constexpr Matrix3(const Vector3& row0, const Vector3& row1, const Vector3& row2)
      : m_{row0.x, row0.y, row0.z, row1.x, row1.y, row1.z, row2.x, row2.y, row2.z} {}

  [[nodiscard]] static constexpr Matrix3 identity() { return {}; }

  /// Rotation by `angle` radians (right-hand rule) about the axis direction
  /// `axis` (need not be unit length), via Rodrigues' formula
  ///   R = I + sin(a) K + (1 - cos(a)) K^2,   K = [axis]_x.
  /// Throws DegenerateGeometryError for a zero axis and InvalidArgumentError for
  /// a non-finite angle.
  [[nodiscard]] static Matrix3 rotationAxisAngle(const Vector3& axis, double angle) {
    if (!std::isfinite(angle)) throw InvalidArgumentError("rotationAxisAngle: non-finite angle");
    const Vector3 u = normalized(axis);
    const double s = std::sin(angle);
    const double c = std::cos(angle);
    const double t = 1.0 - c;
    return Matrix3{{t * u.x * u.x + c, t * u.x * u.y - s * u.z, t * u.x * u.z + s * u.y},
                   {t * u.x * u.y + s * u.z, t * u.y * u.y + c, t * u.y * u.z - s * u.x},
                   {t * u.x * u.z - s * u.y, t * u.y * u.z + s * u.x, t * u.z * u.z + c}};
  }

  [[nodiscard]] constexpr double operator()(int row, int col) const {
    return m_[static_cast<std::size_t>(row * 3 + col)];
  }

  [[nodiscard]] constexpr Vector3 row(int r) const {
    return {(*this)(r, 0), (*this)(r, 1), (*this)(r, 2)};
  }
  [[nodiscard]] constexpr Vector3 col(int c) const {
    return {(*this)(0, c), (*this)(1, c), (*this)(2, c)};
  }

  [[nodiscard]] constexpr Matrix3 transposed() const {
    return Matrix3{col(0), col(1), col(2)};
  }

  [[nodiscard]] constexpr double determinant() const {
    return dot(row(0), cross(row(1), row(2)));
  }

  /// True if R^T R = I (entrywise within `tol`) and det(R) = +1 (within `tol`).
  [[nodiscard]] bool isRotation(double tol = 1e-9) const;

 private:
  std::array<double, 9> m_;
};

[[nodiscard]] constexpr Vector3 operator*(const Matrix3& a, const Vector3& v) {
  return {dot(a.row(0), v), dot(a.row(1), v), dot(a.row(2), v)};
}

[[nodiscard]] constexpr Matrix3 operator*(const Matrix3& a, const Matrix3& b) {
  return Matrix3{{dot(a.row(0), b.col(0)), dot(a.row(0), b.col(1)), dot(a.row(0), b.col(2))},
                 {dot(a.row(1), b.col(0)), dot(a.row(1), b.col(1)), dot(a.row(1), b.col(2))},
                 {dot(a.row(2), b.col(0)), dot(a.row(2), b.col(1)), dot(a.row(2), b.col(2))}};
}

inline bool Matrix3::isRotation(double tol) const {
  for (const double e : m_) {
    if (!std::isfinite(e)) return false;
  }
  const Matrix3 p = transposed() * (*this);
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      if (std::abs(p(r, c) - (r == c ? 1.0 : 0.0)) > tol) return false;
    }
  }
  return std::abs(determinant() - 1.0) <= tol;
}

}  // namespace geoconstraint
