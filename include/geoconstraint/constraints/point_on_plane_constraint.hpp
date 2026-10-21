#pragma once

#include "geoconstraint/constraints/constraint.hpp"
#include "geoconstraint/geometry/plane3.hpp"

namespace geoconstraint {

/// A point lies on a fixed plane.   r = n . (p - p0)   (1 residual)
/// `n` is the plane's unit normal and `p0` its origin, so the residual is the
/// signed distance of the point from the plane (length units).
class PointOnPlaneConstraint final : public Constraint {
 public:
  PointOnPlaneConstraint(PointId p, const Plane3& plane) : p_(p), plane_(plane) {}

  [[nodiscard]] std::string_view typeName() const noexcept override { return "PointOnPlane"; }
  [[nodiscard]] std::size_t residualCount() const noexcept override { return 1; }
  [[nodiscard]] std::vector<PointId> referencedPoints() const override { return {p_}; }
  void evaluate(std::span<const Point3> points, std::span<double> residuals) const override;

 private:
  PointId p_;
  Plane3 plane_;
};

}  // namespace geoconstraint
