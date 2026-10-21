#pragma once

#include "geoconstraint/constraints/constraint.hpp"

namespace geoconstraint {

/// Two line directions are parallel.   r = u x v   (3 residuals)
///
/// u and v are the unit directions of the two `DirectionRef`s. |u x v| = sin of
/// the angle between them, so the constraint is satisfied by BOTH parallel and
/// antiparallel directions (use `AngleConstraint` with theta = 0 to demand the
/// same orientation). Near a satisfied configuration (u = +/-v) the first-order
/// variation of u x v lies in the 2-D plane orthogonal to u, so the 3 x n
/// Jacobian has rank 2 there: one residual row is redundant and the
/// diagnostics report it as such. (Away from parallel the rank is 3.)
/// A zero-length direction raises DegenerateGeometryError - it is never
/// considered parallel to anything.
class ParallelConstraint final : public Constraint {
 public:
  ParallelConstraint(DirectionRef a, DirectionRef b) : a_(a), b_(b) {}

  [[nodiscard]] std::string_view typeName() const noexcept override { return "Parallel"; }
  [[nodiscard]] std::size_t residualCount() const noexcept override { return 3; }
  [[nodiscard]] std::vector<PointId> referencedPoints() const override {
    return {a_.from, a_.to, b_.from, b_.to};
  }
  void evaluate(std::span<const Point3> points, std::span<double> residuals) const override;

 private:
  DirectionRef a_;
  DirectionRef b_;
};

}  // namespace geoconstraint
