#pragma once

#include "geoconstraint/constraints/constraint.hpp"

namespace geoconstraint {

/// Two line directions are perpendicular.   r = u . v   (1 residual)
/// u, v are unit directions, so r = cos(angle between them). A zero-length
/// direction raises DegenerateGeometryError (u . v = 0 would otherwise be
/// satisfied trivially).
class PerpendicularConstraint final : public Constraint {
 public:
  PerpendicularConstraint(DirectionRef a, DirectionRef b) : a_(a), b_(b) {}

  [[nodiscard]] std::string_view typeName() const noexcept override { return "Perpendicular"; }
  [[nodiscard]] std::size_t residualCount() const noexcept override { return 1; }
  [[nodiscard]] std::vector<PointId> referencedPoints() const override {
    return {a_.from, a_.to, b_.from, b_.to};
  }
  void evaluate(std::span<const Point3> points, std::span<double> residuals) const override;

 private:
  DirectionRef a_;
  DirectionRef b_;
};

}  // namespace geoconstraint
