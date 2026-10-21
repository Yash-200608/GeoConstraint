// Tour of the geometry layer. Every printed value is also checked against the
// analytically known answer; the program exits non-zero if any check fails.

#include <cmath>
#include <iomanip>
#include <iostream>

#include "geoconstraint/geometry/distances.hpp"
#include "geoconstraint/geometry/intersections.hpp"
#include "geoconstraint/geometry/projections.hpp"
#include "geoconstraint/math/rigid_transform3.hpp"

using namespace geoconstraint;

namespace {

int failures = 0;

void check(const char* what, bool condition) {
  std::cout << (condition ? "  [ok]   " : "  [FAIL] ") << what << "\n";
  if (!condition) ++failures;
}

void print(const char* name, const Point3& p) {
  std::cout << "  " << name << " = (" << p.x << ", " << p.y << ", " << p.z << ")\n";
}

}  // namespace

int main() {
  std::cout << std::fixed << std::setprecision(6);

  std::cout << "Vectors\n";
  const Vector3 u{1, 0, 0};
  const Vector3 v{0, 2, 0};
  const Vector3 w = cross(u, v);
  std::cout << "  u x v = (" << w.x << ", " << w.y << ", " << w.z << ")\n";
  check("u x v = (0,0,2)", approxEqual(w, Vector3{0, 0, 2}));
  check("u perpendicular to v", isPerpendicular(u, v));
  check("angle(u, v) = pi/2", std::abs(angleBetween(u, v) - std::acos(-1.0) / 2) < 1e-12);

  std::cout << "Projection onto an infinite line / plane\n";
  const Line3 xAxis = Line3::fromPointDirection({0, 0, 0}, {1, 0, 0});
  const Plane3 xy = Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1});
  const Point3 p{3, 4, 5};
  const Point3 onLine = projectPointOnLine(p, xAxis);
  const Point3 onPlane = projectPointOnPlane(p, xy);
  print("projection on x-axis", onLine);
  print("projection on z=0   ", onPlane);
  check("line projection = (3,0,0)", approxEqual(onLine, Point3{3, 0, 0}));
  check("plane projection = (3,4,0)", approxEqual(onPlane, Point3{3, 4, 0}));
  std::cout << "  distance to x-axis = " << distancePointToLine(p, xAxis)
            << ", distance to plane = " << distancePointToPlane(p, xy) << "\n";
  check("distance to line = sqrt(41)", approxEqual(distancePointToLine(p, xAxis), std::sqrt(41.0)));

  std::cout << "Line-plane intersection\n";
  const Line3 slanted = Line3::fromTwoPoints({0, 0, 2}, {1, 1, 1});
  const auto hit = intersectLinePlane(slanted, xy);
  check("slanted line crosses z=0 once", hit.relation == LinePlaneRelation::UniquePoint);
  if (hit.point) print("hit", *hit.point);
  check("hit = (2,2,0)", hit.point && approxEqual(*hit.point, Point3{2, 2, 0}));
  const Line3 parallelAbove = Line3::fromPointDirection({0, 0, 1}, {1, 0, 0});
  check("parallel line above plane: no intersection",
        intersectLinePlane(parallelAbove, xy).relation == LinePlaneRelation::NoIntersection);
  check("line inside plane detected",
        intersectLinePlane(Line3::fromPointDirection({0, 0, 0}, {1, 1, 0}), xy).relation ==
            LinePlaneRelation::LineInPlane);

  std::cout << "Closest points of two skew lines\n";
  const Line3 lineA = Line3::fromPointDirection({0, 0, 0}, {1, 0, 0});
  const Line3 lineB = Line3::fromPointDirection({0, 0, 3}, {0, 1, 0});
  const auto closest = closestPointsBetweenLines(lineA, lineB);
  print("on A", closest.pointOnA);
  print("on B", closest.pointOnB);
  std::cout << "  distance = " << closest.distance << "\n";
  check("lines are skew at distance 3",
        closest.relation == LineLineRelation::Skew && approxEqual(closest.distance, 3.0));

  std::cout << "Rigid transforms\n";
  const double halfPi = std::acos(-1.0) / 2;
  const RigidTransform3 rot = RigidTransform3::rotation({0, 0, 1}, halfPi);
  const RigidTransform3 shift = RigidTransform3::translation({10, 0, 0});
  const RigidTransform3 both = shift * rot;  // rotate first, then translate
  const Point3 moved = both.apply(Point3{1, 0, 0});
  print("(1,0,0) rotated 90 deg about z, then moved +10 in x", moved);
  check("result = (10,1,0)", approxEqual(moved, Point3{10, 1, 0}));
  check("inverse restores the point",
        approxEqual(both.inverse().apply(moved), Point3{1, 0, 0}));

  std::cout << (failures == 0 ? "All checks passed\n" : "Some checks FAILED\n");
  return failures == 0 ? 0 : 1;
}
