#pragma once

#include "geoconstraint/constraints/constraint.hpp"

namespace geoconstraint {

/// Two points coincide.   r = p - q   (3 residuals: x, y, z)
class CoincidenceConstraint final : public Constraint {
 public:
  CoincidenceConstraint(PointId p, PointId q) : p_(p), q_(q) {}

  [[nodiscard]] std::string_view typeName() const noexcept override { return "Coincident"; }
  [[nodiscard]] std::size_t residualCount() const noexcept override { return 3; }
  [[nodiscard]] std::vector<PointId> referencedPoints() const override { return {p_, q_}; }
  void evaluate(std::span<const Point3> points, std::span<double> residuals) const override;

 private:
  PointId p_;
  PointId q_;
};

}  // namespace geoconstraint
