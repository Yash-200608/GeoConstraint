#pragma once

#include <optional>

#include "geoconstraint/geometry/line3.hpp"
#include "geoconstraint/geometry/plane3.hpp"

namespace geoconstraint {

enum class LinePlaneRelation {
  UniquePoint,     ///< the line crosses the plane at exactly one point
  NoIntersection,  ///< parallel to the plane and not contained in it
  LineInPlane,     ///< parallel and contained: every line point is common
  Invalid          ///< non-finite input or non-finite result
};

struct LinePlaneIntersection {
  LinePlaneRelation relation{LinePlaneRelation::Invalid};
  /// Present only for UniquePoint.
  std::optional<Point3> point;
  /// Line parameter t of `point` (point = origin + t * direction); present only
  /// for UniquePoint.
  std::optional<double> parameter;
};

/// Intersection of an INFINITE line with an INFINITE plane.
///
/// With unit direction d and unit normal n, denom = n . d equals the sine of the
/// angle between the line and the plane. If |denom| <= tol.angular the line is
/// treated as parallel: it lies in the plane when the signed distance of its
/// origin is within the scale-aware length tolerance, otherwise there is no
/// intersection. Otherwise t = n . (planeOrigin - lineOrigin) / denom.
///
/// A line whose tilt is just above the threshold is a legitimate crossing but
/// may produce a very distant point; see docs/numerical_robustness.md.
[[nodiscard]] LinePlaneIntersection intersectLinePlane(
    const Line3& line, const Plane3& plane, const Tolerances& tol = kDefaultTolerances);

}  // namespace geoconstraint
