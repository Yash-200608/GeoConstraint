#include <algorithm>
#include <cmath>

#include "geoconstraint/geometry/intersections.hpp"

namespace geoconstraint {

LinePlaneIntersection intersectLinePlane(const Line3& line, const Plane3& plane,
                                         const Tolerances& tol) {
  const Vector3& d = line.direction();
  const Vector3& n = plane.normal();
  const double denom = dot(n, d);
  const double offset = plane.signedDistance(line.origin());  // n . (o_line - o_plane)
  if (!std::isfinite(denom) || !std::isfinite(offset)) return {};

  if (std::abs(denom) <= tol.angular) {
    const double scale = std::max(norm(line.origin().positionVector()),
                                  norm(plane.origin().positionVector()));
    const bool inPlane = std::abs(offset) <= lengthTolerance(tol, scale);
    return {inPlane ? LinePlaneRelation::LineInPlane : LinePlaneRelation::NoIntersection,
            std::nullopt, std::nullopt};
  }

  const double t = -offset / denom;
  const Point3 hit = line.pointAt(t);
  if (!std::isfinite(t) || !hit.isFinite()) return {};
  return {LinePlaneRelation::UniquePoint, hit, t};
}

}  // namespace geoconstraint
