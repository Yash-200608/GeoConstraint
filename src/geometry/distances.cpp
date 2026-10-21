#include "geoconstraint/geometry/distances.hpp"

#include <algorithm>
#include <cmath>

namespace geoconstraint {

double distancePointToLine(const Point3& p, const Line3& line) {
  return norm(cross(p - line.origin(), line.direction()));
}

double distancePointToPlane(const Point3& p, const Plane3& plane) {
  return std::abs(plane.signedDistance(p));
}

LineLineClosest closestPointsBetweenLines(const Line3& a, const Line3& b, const Tolerances& tol) {
  LineLineClosest out;
  const Vector3 w = b.origin() - a.origin();
  const Vector3 n = cross(a.direction(), b.direction());
  const double nn = squaredNorm(n);
  const double sinAngle = std::sqrt(nn);  // both directions are unit vectors
  const double scale = std::max(norm(a.origin().positionVector()),
                                norm(b.origin().positionVector()));
  const double lengthTol = lengthTolerance(tol, scale);

  if (sinAngle <= tol.angular) {
    // Parallel: the distance is that of B's origin to line A; the closest-point
    // pair is not unique, so report the pair (A origin, its projection on B).
    out.parameterA = 0.0;
    out.pointOnA = a.origin();
    out.parameterB = b.parameterOf(a.origin());
    out.pointOnB = b.pointAt(out.parameterB);
    out.distance = distance(out.pointOnA, out.pointOnB);
    out.relation = out.distance <= lengthTol ? LineLineRelation::Coincident
                                             : LineLineRelation::Parallel;
  } else {
    out.parameterA = dot(cross(w, b.direction()), n) / nn;
    out.parameterB = dot(cross(w, a.direction()), n) / nn;
    out.pointOnA = a.pointAt(out.parameterA);
    out.pointOnB = b.pointAt(out.parameterB);
    out.distance = std::abs(dot(w, n)) / sinAngle;
    out.relation = out.distance <= lengthTol ? LineLineRelation::Intersecting
                                             : LineLineRelation::Skew;
  }

  if (!out.pointOnA.isFinite() || !out.pointOnB.isFinite() || !std::isfinite(out.distance) ||
      !std::isfinite(out.parameterA) || !std::isfinite(out.parameterB)) {
    return LineLineClosest{};  // relation == Invalid
  }
  return out;
}

}  // namespace geoconstraint
