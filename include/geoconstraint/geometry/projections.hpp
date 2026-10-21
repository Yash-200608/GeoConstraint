#pragma once

#include "geoconstraint/geometry/line3.hpp"
#include "geoconstraint/geometry/plane3.hpp"

namespace geoconstraint {

/// Orthogonal projection of `p` onto the INFINITE line:
///   o + ((p - o) . d) d.
/// The result may lie anywhere on the unbounded line (no clamping to a segment).
[[nodiscard]] inline Point3 projectPointOnLine(const Point3& p, const Line3& line) {
  return line.pointAt(line.parameterOf(p));
}

/// Orthogonal projection of `p` onto the INFINITE plane:
///   p - (n . (p - o)) n.
[[nodiscard]] inline Point3 projectPointOnPlane(const Point3& p, const Plane3& plane) {
  return p - plane.normal() * plane.signedDistance(p);
}

}  // namespace geoconstraint
