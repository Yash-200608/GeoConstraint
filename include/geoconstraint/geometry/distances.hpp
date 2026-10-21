#pragma once

#include "geoconstraint/geometry/line3.hpp"
#include "geoconstraint/geometry/plane3.hpp"

namespace geoconstraint {

// Point-point distance is `distance(Point3, Point3)` in point3.hpp.

/// Distance from `p` to the INFINITE line: |(p - o) x d| (d is unit length).
[[nodiscard]] double distancePointToLine(const Point3& p, const Line3& line);

/// Unsigned distance from `p` to the INFINITE plane: |n . (p - o)|.
[[nodiscard]] double distancePointToPlane(const Point3& p, const Plane3& plane);

/// Signed distance; positive on the side the plane normal points to.
[[nodiscard]] inline double signedDistancePointToPlane(const Point3& p, const Plane3& plane) {
  return plane.signedDistance(p);
}

/// Relationship of two infinite lines. Classification uses
/// `Tolerances::angular` (on sin of the angle between the unit directions) for
/// "parallel" and the scale-aware length tolerance for "touching".
enum class LineLineRelation {
  Skew,          ///< not parallel, not touching: unique common perpendicular
  Intersecting,  ///< not parallel, closest distance within length tolerance
  Parallel,      ///< parallel, distinct (closest points are not unique)
  Coincident,    ///< the same line (within tolerance)
  Invalid        ///< non-finite result (overflow); no geometry is returned
};

struct LineLineClosest {
  LineLineRelation relation{LineLineRelation::Invalid};
  /// Closest points and their line parameters. Unique for Skew/Intersecting.
  /// For Parallel / Coincident they are a *representative* pair (point on A at
  /// A's origin and its projection on B); the distance is still exact.
  Point3 pointOnA{};
  Point3 pointOnB{};
  double parameterA{0.0};
  double parameterB{0.0};
  /// Minimum distance between the two infinite lines (0 if intersecting or
  /// coincident up to tolerance); meaningless for Invalid.
  double distance{0.0};
};

/// Closest points / distance between two INFINITE lines. For non-parallel
/// lines, with n = dA x dB and w = oB - oA:
///   tA = ((w x dB) . n) / |n|^2,  tB = ((w x dA) . n) / |n|^2,
///   distance = |w . n| / |n|.
/// Nearly parallel lines (sin(angle) > tol.angular but tiny) are classified as
/// non-parallel and yield large parameters; the distance remains accurate to
/// roughly eps / sin(angle) relative error. See docs/numerical_robustness.md.
[[nodiscard]] LineLineClosest closestPointsBetweenLines(
    const Line3& a, const Line3& b, const Tolerances& tol = kDefaultTolerances);

}  // namespace geoconstraint
