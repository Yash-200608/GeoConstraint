#include <cmath>
#include <numbers>
#include <string>

#include "geoconstraint/constraints/angle_constraint.hpp"
#include "geoconstraint/constraints/coincidence_constraint.hpp"
#include "geoconstraint/constraints/constraint.hpp"
#include "geoconstraint/constraints/distance_constraint.hpp"
#include "geoconstraint/constraints/fixed_point_constraint.hpp"
#include "geoconstraint/constraints/parallel_constraint.hpp"
#include "geoconstraint/constraints/perpendicular_constraint.hpp"
#include "geoconstraint/constraints/point_on_plane_constraint.hpp"

namespace geoconstraint {

// ---------------------------------------------------------------------------
// Constraint base helpers
// ---------------------------------------------------------------------------

std::vector<double> Constraint::residuals(std::span<const Point3> points) const {
  std::vector<double> out(residualCount());
  evaluate(points, out);
  return out;
}

const Point3& Constraint::pointAt(std::span<const Point3> points, PointId id) {
  if (id.index >= points.size()) {
    throw InvalidArgumentError("Constraint: point id " + std::to_string(id.index) +
                               " out of range (" + std::to_string(points.size()) + " points)");
  }
  return points[id.index];
}

Vector3 Constraint::unitDirection(std::span<const Point3> points, DirectionRef ref) {
  const Vector3 d = pointAt(points, ref.to) - pointAt(points, ref.from);
  if (auto u = tryNormalized(d)) return *u;
  throw DegenerateGeometryError("Constraint: direction from point " + std::to_string(ref.from.index) +
                                " to point " + std::to_string(ref.to.index) +
                                " is zero-length or non-finite");
}

void Constraint::requireSize(std::span<double> residuals, std::size_t expected) {
  if (residuals.size() != expected) {
    throw InvalidArgumentError("Constraint: residual buffer has size " +
                               std::to_string(residuals.size()) + ", expected " +
                               std::to_string(expected));
  }
}

// ---------------------------------------------------------------------------
// Coincident: r = p - q
// ---------------------------------------------------------------------------

void CoincidenceConstraint::evaluate(std::span<const Point3> points,
                                     std::span<double> residuals) const {
  requireSize(residuals, 3);
  const Vector3 r = pointAt(points, p_) - pointAt(points, q_);
  residuals[0] = r.x;
  residuals[1] = r.y;
  residuals[2] = r.z;
}

// ---------------------------------------------------------------------------
// FixedPoint: r = p - p0
// ---------------------------------------------------------------------------

FixedPointConstraint::FixedPointConstraint(PointId p, const Point3& target)
    : p_(p), target_(target) {
  if (!target.isFinite()) throw InvalidArgumentError("FixedPointConstraint: non-finite target");
}

void FixedPointConstraint::evaluate(std::span<const Point3> points,
                                    std::span<double> residuals) const {
  requireSize(residuals, 3);
  const Vector3 r = pointAt(points, p_) - target_;
  residuals[0] = r.x;
  residuals[1] = r.y;
  residuals[2] = r.z;
}

// ---------------------------------------------------------------------------
// Distance: r = |p - q| - d
// ---------------------------------------------------------------------------

DistanceConstraint::DistanceConstraint(PointId p, PointId q, double distance)
    : p_(p), q_(q), distance_(distance) {
  if (!std::isfinite(distance) || distance < 0.0) {
    throw InvalidArgumentError("DistanceConstraint: distance must be finite and >= 0");
  }
}

void DistanceConstraint::evaluate(std::span<const Point3> points,
                                  std::span<double> residuals) const {
  requireSize(residuals, 1);
  residuals[0] = distance(pointAt(points, p_), pointAt(points, q_)) - distance_;
}

// ---------------------------------------------------------------------------
// PointOnPlane: r = n . (p - p0)
// ---------------------------------------------------------------------------

void PointOnPlaneConstraint::evaluate(std::span<const Point3> points,
                                      std::span<double> residuals) const {
  requireSize(residuals, 1);
  residuals[0] = plane_.signedDistance(pointAt(points, p_));
}

// ---------------------------------------------------------------------------
// Parallel: r = u x v
// ---------------------------------------------------------------------------

void ParallelConstraint::evaluate(std::span<const Point3> points,
                                  std::span<double> residuals) const {
  requireSize(residuals, 3);
  const Vector3 r = cross(unitDirection(points, a_), unitDirection(points, b_));
  residuals[0] = r.x;
  residuals[1] = r.y;
  residuals[2] = r.z;
}

// ---------------------------------------------------------------------------
// Perpendicular: r = u . v
// ---------------------------------------------------------------------------

void PerpendicularConstraint::evaluate(std::span<const Point3> points,
                                       std::span<double> residuals) const {
  requireSize(residuals, 1);
  residuals[0] = dot(unitDirection(points, a_), unitDirection(points, b_));
}

// ---------------------------------------------------------------------------
// Angle: r = u . v - cos(theta)
// ---------------------------------------------------------------------------

AngleConstraint::AngleConstraint(DirectionRef a, DirectionRef b, double angleRadians)
    : a_(a), b_(b), cosTarget_(0.0) {
  if (!std::isfinite(angleRadians) || angleRadians < 0.0 || angleRadians > std::numbers::pi) {
    throw InvalidArgumentError("AngleConstraint: angle must be finite and in [0, pi]");
  }
  cosTarget_ = std::cos(angleRadians);
}

void AngleConstraint::evaluate(std::span<const Point3> points, std::span<double> residuals) const {
  requireSize(residuals, 1);
  residuals[0] = dot(unitDirection(points, a_), unitDirection(points, b_)) - cosTarget_;
}

}  // namespace geoconstraint
