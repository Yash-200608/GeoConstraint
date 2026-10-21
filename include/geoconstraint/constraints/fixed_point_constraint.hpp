#pragma once

#include "geoconstraint/constraints/constraint.hpp"

namespace geoconstraint {

/// A point is pinned to a target location.   r = p - p0   (3 residuals)
///
/// This is a *soft* constraint like every other: the solver drives it to zero
/// within its tolerance rather than removing the point from the unknowns.
class FixedPointConstraint final : public Constraint {
 public:
  /// Throws InvalidArgumentError if `target` is non-finite.
  FixedPointConstraint(PointId p, const Point3& target);

  [[nodiscard]] std::string_view typeName() const noexcept override { return "FixedPoint"; }
  [[nodiscard]] std::size_t residualCount() const noexcept override { return 3; }
  [[nodiscard]] std::vector<PointId> referencedPoints() const override { return {p_}; }
  void evaluate(std::span<const Point3> points, std::span<double> residuals) const override;

 private:
  PointId p_;
  Point3 target_;
};

}  // namespace geoconstraint
