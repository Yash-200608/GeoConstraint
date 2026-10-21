#pragma once

#include "geoconstraint/constraints/constraint.hpp"

namespace geoconstraint {

/// Fixed angle between two DIRECTED line directions.   r = u . v - cos(theta)
/// (1 residual).
///
/// Convention
///  - u = unit(P[a.to] - P[a.from]), v = unit(P[b.to] - P[b.from]); the order of
///    the two points in each `DirectionRef` matters (reversing one direction
///    changes theta to pi - theta).
///  - theta is the *unsigned* angle between u and v, in [0, pi]. There is no
///    orientation/handedness: the angle is the same for (a, b) and (b, a).
///  - theta = 0 selects parallel same-direction, theta = pi antiparallel,
///    theta = pi/2 perpendicular. Antiparallel directions are therefore NOT
///    accepted by a theta = 0 constraint (u . v = -1 != 1).
///  - Because the residual is cos-based, its gradient vanishes at theta = 0 and
///    theta = pi (d cos / d theta = -sin theta = 0). The system is then
///    first-order singular there: the Jacobian loses rank and convergence from
///    a perturbed start is slow (linear, not quadratic). See docs/limitations.md.
///  - A zero-length direction raises DegenerateGeometryError.
class AngleConstraint final : public Constraint {
 public:
  /// Throws InvalidArgumentError unless `angleRadians` is finite and in [0, pi].
  AngleConstraint(DirectionRef a, DirectionRef b, double angleRadians);

  [[nodiscard]] std::string_view typeName() const noexcept override { return "Angle"; }
  [[nodiscard]] std::size_t residualCount() const noexcept override { return 1; }
  [[nodiscard]] std::vector<PointId> referencedPoints() const override {
    return {a_.from, a_.to, b_.from, b_.to};
  }
  void evaluate(std::span<const Point3> points, std::span<double> residuals) const override;

 private:
  DirectionRef a_;
  DirectionRef b_;
  double cosTarget_;
};

}  // namespace geoconstraint
