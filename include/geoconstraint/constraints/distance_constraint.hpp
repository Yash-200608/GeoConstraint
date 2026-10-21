#pragma once

#include "geoconstraint/constraints/constraint.hpp"

namespace geoconstraint {

/// Fixed distance between two points.   r = |p - q| - d   (1 residual)
///
/// Note: |p - q| is not differentiable at p == q. A finite-difference Jacobian
/// there is meaningless (central differences give a zero gradient), so avoid
/// starting the solver with the two points exactly coincident when d > 0.
class DistanceConstraint final : public Constraint {
 public:
  /// Throws InvalidArgumentError unless `distance` is finite and >= 0.
  DistanceConstraint(PointId p, PointId q, double distance);

  [[nodiscard]] std::string_view typeName() const noexcept override { return "Distance"; }
  [[nodiscard]] std::size_t residualCount() const noexcept override { return 1; }
  [[nodiscard]] std::vector<PointId> referencedPoints() const override { return {p_, q_}; }
  void evaluate(std::span<const Point3> points, std::span<double> residuals) const override;

 private:
  PointId p_;
  PointId q_;
  double distance_;
};

}  // namespace geoconstraint
