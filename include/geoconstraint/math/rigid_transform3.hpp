#pragma once

#include "geoconstraint/math/matrix3.hpp"
#include "geoconstraint/math/point3.hpp"

namespace geoconstraint {

/// Proper rigid motion x -> R x + t with R a rotation (R^T R = I, det R = +1).
/// Reflections and scalings are rejected at construction. Points are rotated
/// and translated; vectors are only rotated.
class RigidTransform3 {
 public:
  /// Identity transform.
  RigidTransform3() = default;

  /// Throws InvalidArgumentError if `rotation` is not a proper rotation
  /// (tolerance 1e-9) or `translation` is non-finite.
  RigidTransform3(const Matrix3& rotation, const Vector3& translation)
      : rotation_(rotation), translation_(translation) {
    if (!rotation.isRotation()) {
      throw InvalidArgumentError("RigidTransform3: matrix is not a proper rotation");
    }
    if (!translation.isFinite()) {
      throw InvalidArgumentError("RigidTransform3: non-finite translation");
    }
  }

  [[nodiscard]] static RigidTransform3 translation(const Vector3& t) {
    return RigidTransform3(Matrix3::identity(), t);
  }
  /// Rotation about an axis through the origin.
  [[nodiscard]] static RigidTransform3 rotation(const Vector3& axis, double angle) {
    return RigidTransform3(Matrix3::rotationAxisAngle(axis, angle), Vector3{});
  }
  /// Rotation about an axis through `center` (the center is a fixed point).
  [[nodiscard]] static RigidTransform3 rotationAbout(const Point3& center, const Vector3& axis,
                                                     double angle) {
    const Matrix3 r = Matrix3::rotationAxisAngle(axis, angle);
    return RigidTransform3(r, center.positionVector() - r * center.positionVector());
  }

  [[nodiscard]] const Matrix3& rotationMatrix() const noexcept { return rotation_; }
  [[nodiscard]] const Vector3& translationVector() const noexcept { return translation_; }

  [[nodiscard]] Point3 apply(const Point3& p) const {
    const Vector3 v = rotation_ * p.positionVector() + translation_;
    return {v.x, v.y, v.z};
  }
  /// Vectors are directions: translation does not apply.
  [[nodiscard]] Vector3 apply(const Vector3& v) const { return rotation_ * v; }

  /// (a * b).apply(p) == a.apply(b.apply(p)): `b` acts first.
  [[nodiscard]] RigidTransform3 operator*(const RigidTransform3& b) const {
    return RigidTransform3(rotation_ * b.rotation_, rotation_ * b.translation_ + translation_);
  }

  /// Inverse: x -> R^T (x - t).
  [[nodiscard]] RigidTransform3 inverse() const {
    const Matrix3 rt = rotation_.transposed();
    return RigidTransform3(rt, -(rt * translation_));
  }

 private:
  Matrix3 rotation_{};
  Vector3 translation_{};
};

}  // namespace geoconstraint
