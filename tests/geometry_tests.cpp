#include <cmath>

#include "geoconstraint/geometry/distances.hpp"
#include "geoconstraint/geometry/intersections.hpp"
#include "geoconstraint/geometry/line3.hpp"
#include "geoconstraint/geometry/plane3.hpp"
#include "geoconstraint/geometry/projections.hpp"
#include "geoconstraint/math/rigid_transform3.hpp"
#include "test_utils.hpp"

using namespace geoconstraint;

TEST_CASE("Line3 and Plane3 construction validates input", "[geometry][errors]") {
  REQUIRE_THROWS_AS(Line3::fromPointDirection({0, 0, 0}, {0, 0, 0}), DegenerateGeometryError);
  REQUIRE_THROWS_AS(Line3::fromTwoPoints({1, 2, 3}, {1, 2, 3}), DegenerateGeometryError);
  REQUIRE_THROWS_AS(Line3::fromPointDirection({std::nan(""), 0, 0}, {1, 0, 0}),
                    InvalidArgumentError);
  REQUIRE_THROWS_AS(Plane3::fromPointNormal({0, 0, 0}, {0, 0, 0}), DegenerateGeometryError);
  REQUIRE_THROWS_AS(Plane3::fromPointNormal({0, 0, 0}, {std::nan(""), 0, 1}),
                    DegenerateGeometryError);
  REQUIRE_THROWS_AS(Plane3::fromPointNormal({0, 0, std::nan("")}, {0, 0, 1}),
                    InvalidArgumentError);
  // three points: collinear and coincident
  REQUIRE_THROWS_AS(Plane3::fromThreePoints({0, 0, 0}, {1, 1, 1}, {2, 2, 2}),
                    DegenerateGeometryError);
  REQUIRE_THROWS_AS(Plane3::fromThreePoints({0, 0, 0}, {0, 0, 0}, {1, 2, 3}),
                    DegenerateGeometryError);
}

TEST_CASE("Line3 and Plane3 normalize their direction/normal", "[geometry]") {
  const Line3 line = Line3::fromPointDirection({1, 1, 1}, {0, 0, 5});
  REQUIRE_NEAR(norm(line.direction()), 1.0, 1e-15);
  REQUIRE(line.pointAt(2.0) == Point3{1, 1, 3});
  const Plane3 plane = Plane3::fromPointNormal({0, 0, 4}, {0, 0, -9});
  REQUIRE_NEAR(norm(plane.normal()), 1.0, 1e-15);
  REQUIRE_NEAR(plane.signedDistance({7, 7, 6}), -2.0, 1e-15);  // normal points to -z
}

TEST_CASE("Plane from three points follows the right-hand rule", "[geometry]") {
  const Plane3 p = Plane3::fromThreePoints({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
  REQUIRE_POINT_NEAR(p.normal(), (Vector3{0, 0, 1}), 1e-15);
  const Plane3 q = Plane3::fromThreePoints({0, 0, 0}, {0, 1, 0}, {1, 0, 0});
  REQUIRE_POINT_NEAR(q.normal(), (Vector3{0, 0, -1}), 1e-15);
}

TEST_CASE("Projection onto an infinite line", "[geometry][projection]") {
  const Line3 xAxis = Line3::fromPointDirection({0, 0, 0}, {1, 0, 0});
  REQUIRE_POINT_NEAR(projectPointOnLine({3, 4, 5}, xAxis), (Point3{3, 0, 0}), 1e-15);
  // The projection may lie beyond any "segment" - lines are infinite.
  REQUIRE_POINT_NEAR(projectPointOnLine({-1000, 4, 5}, xAxis), (Point3{-1000, 0, 0}), 1e-12);

  const Line3 diag = Line3::fromPointDirection({1, 1, 0}, {1, 1, 0});
  REQUIRE_POINT_NEAR(projectPointOnLine({3, 1, 7}, diag), (Point3{2, 2, 0}), 1e-14);

  // The residual p - proj(p) is perpendicular to the line.
  const Point3 p{2.5, -1.5, 4.0};
  const Line3 oblique = Line3::fromPointDirection({1, 2, 3}, {1, -2, 0.5});
  const Vector3 residual = p - projectPointOnLine(p, oblique);
  REQUIRE_NEAR(dot(residual, oblique.direction()), 0.0, 1e-14);
  // A point already on the line is its own projection.
  REQUIRE_POINT_NEAR(projectPointOnLine(oblique.pointAt(-3.7), oblique), oblique.pointAt(-3.7),
                     1e-14);
}

TEST_CASE("Projection onto an infinite plane", "[geometry][projection]") {
  const Plane3 xy = Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1});
  REQUIRE_POINT_NEAR(projectPointOnPlane({3, 4, 5}, xy), (Point3{3, 4, 0}), 1e-15);
  const Plane3 shifted = Plane3::fromPointNormal({0, 0, 2}, {0, 0, 1});
  REQUIRE_POINT_NEAR(projectPointOnPlane({3, 4, -5}, shifted), (Point3{3, 4, 2}), 1e-15);

  const Plane3 tilted = Plane3::fromPointNormal({1, 0, 0}, {1, 1, 0});  // x + y = 1
  const Point3 proj = projectPointOnPlane({3, 3, 9}, tilted);
  REQUIRE_POINT_NEAR(proj, (Point3{0.5, 0.5, 9}), 1e-14);
  REQUIRE_NEAR(tilted.signedDistance(proj), 0.0, 1e-14);
  // Idempotent
  REQUIRE_POINT_NEAR(projectPointOnPlane(proj, tilted), proj, 1e-14);
}

TEST_CASE("Point-to-line and point-to-plane distances", "[geometry][distance]") {
  const Line3 xAxis = Line3::fromPointDirection({0, 0, 0}, {1, 0, 0});
  REQUIRE_NEAR(distancePointToLine({3, 4, 0}, xAxis), 4.0, 1e-15);
  REQUIRE_NEAR(distancePointToLine({3, 3, 4}, xAxis), 5.0, 1e-15);
  REQUIRE_NEAR(distancePointToLine({9, 0, 0}, xAxis), 0.0, 1e-15);
  // Line through (0,0,0) with direction (1,1,0): distance of (1,0,0) is 1/sqrt(2)
  const Line3 diag = Line3::fromPointDirection({0, 0, 0}, {1, 1, 0});
  REQUIRE_NEAR(distancePointToLine({1, 0, 0}, diag), std::sqrt(0.5), 1e-15);

  const Plane3 plane = Plane3::fromPointNormal({0, 0, 1}, {0, 0, 2});
  REQUIRE_NEAR(distancePointToPlane({5, 5, 4}, plane), 3.0, 1e-15);
  REQUIRE_NEAR(distancePointToPlane({5, 5, -2}, plane), 3.0, 1e-15);
  REQUIRE_NEAR(signedDistancePointToPlane({5, 5, 4}, plane), 3.0, 1e-15);
  REQUIRE_NEAR(signedDistancePointToPlane({5, 5, -2}, plane), -3.0, 1e-15);
  REQUIRE_NEAR(distance(Point3{1, 2, 3}, Point3{4, 6, 3}), 5.0, 1e-15);
}

TEST_CASE("Line-plane intersection: unique point", "[geometry][intersection]") {
  const Plane3 xy = Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1});
  const auto hit = intersectLinePlane(Line3::fromTwoPoints({0, 0, 2}, {1, 1, 1}), xy);
  REQUIRE(hit.relation == LinePlaneRelation::UniquePoint);
  REQUIRE(hit.point.has_value());
  REQUIRE_POINT_NEAR(*hit.point, (Point3{2, 2, 0}), 1e-14);
  // The reported parameter reproduces the point.
  const Line3 line = Line3::fromTwoPoints({0, 0, 2}, {1, 1, 1});
  REQUIRE_POINT_NEAR(line.pointAt(*hit.parameter), *hit.point, 1e-14);

  // Tilted plane x + y + z = 3, line along (1,2,3)-direction from the origin.
  const Plane3 tilted = Plane3::fromPointNormal({1, 1, 1}, {1, 1, 1});
  const auto h2 = intersectLinePlane(Line3::fromPointDirection({0, 0, 0}, {1, 2, 3}), tilted);
  REQUIRE(h2.relation == LinePlaneRelation::UniquePoint);
  REQUIRE_POINT_NEAR(*h2.point, (Point3{0.5, 1.0, 1.5}), 1e-14);  // t(1+2+3)=3 -> t=0.5
  // A line crossing "behind" its origin still intersects (infinite line).
  const auto h3 = intersectLinePlane(Line3::fromPointDirection({0, 0, 5}, {0, 0, 1}), xy);
  REQUIRE(h3.relation == LinePlaneRelation::UniquePoint);
  REQUIRE_POINT_NEAR(*h3.point, (Point3{0, 0, 0}), 1e-15);
  REQUIRE(*h3.parameter < 0.0);
}

TEST_CASE("Line-plane intersection: parallel and contained lines", "[geometry][intersection]") {
  const Plane3 xy = Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1});
  const auto above = intersectLinePlane(Line3::fromPointDirection({0, 0, 1}, {1, 2, 0}), xy);
  REQUIRE(above.relation == LinePlaneRelation::NoIntersection);
  REQUIRE_FALSE(above.point.has_value());
  REQUIRE_FALSE(above.parameter.has_value());

  const auto inside = intersectLinePlane(Line3::fromPointDirection({5, 5, 0}, {1, 2, 0}), xy);
  REQUIRE(inside.relation == LinePlaneRelation::LineInPlane);
  REQUIRE_FALSE(inside.point.has_value());  // not unique: no single point is reported
}

TEST_CASE("Line-plane intersection: near-parallel configurations", "[geometry][intersection]") {
  const Plane3 xy = Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1});
  // Tilt of 1e-6 rad is far above the 1e-9 threshold: a real, distant crossing.
  const double eps = 1e-6;
  const Line3 shallow =
      Line3::fromPointDirection({0, 0, 1}, {std::cos(eps), 0, -std::sin(eps)});
  const auto hit = intersectLinePlane(shallow, xy);
  REQUIRE(hit.relation == LinePlaneRelation::UniquePoint);
  REQUIRE_THAT(hit.point->x, Catch::Matchers::WithinRel(1.0 / std::tan(eps), 1e-9));
  REQUIRE_NEAR(hit.point->z, 0.0, 1e-6);

  // Tilt of 1e-12 rad is below the threshold: classified as parallel.
  const double tiny = 1e-12;
  const Line3 almostParallelAbove =
      Line3::fromPointDirection({0, 0, 1}, {std::cos(tiny), 0, std::sin(tiny)});
  REQUIRE(intersectLinePlane(almostParallelAbove, xy).relation ==
          LinePlaneRelation::NoIntersection);
  const Line3 almostParallelInside =
      Line3::fromPointDirection({0, 0, 0}, {std::cos(tiny), 0, std::sin(tiny)});
  REQUIRE(intersectLinePlane(almostParallelInside, xy).relation == LinePlaneRelation::LineInPlane);
  // A line 1e-6 above the plane is not "in" the plane.
  const Line3 closeAbove = Line3::fromPointDirection({0, 0, 1e-6}, {1, 0, 0});
  REQUIRE(intersectLinePlane(closeAbove, xy).relation == LinePlaneRelation::NoIntersection);
}

TEST_CASE("Closest points: skew lines", "[geometry][lines]") {
  const Line3 a = Line3::fromPointDirection({0, 0, 0}, {1, 0, 0});
  const Line3 b = Line3::fromPointDirection({1, 2, 3}, {0, 1, 0});
  const auto r = closestPointsBetweenLines(a, b);
  REQUIRE(r.relation == LineLineRelation::Skew);
  REQUIRE_NEAR(r.distance, 3.0, 1e-14);
  REQUIRE_POINT_NEAR(r.pointOnA, (Point3{1, 0, 0}), 1e-14);
  REQUIRE_POINT_NEAR(r.pointOnB, (Point3{1, 0, 3}), 1e-14);
}

TEST_CASE("Closest points of skew lines are mutually perpendicular", "[geometry][lines]") {
  // Independent property check on an oblique configuration (no closed form used).
  const Line3 a = Line3::fromPointDirection({1, -2, 0.5}, {1, 2, 3});
  const Line3 b = Line3::fromPointDirection({-4, 3, 2}, {-2, 1, 1});
  const auto r = closestPointsBetweenLines(a, b);
  REQUIRE(r.relation == LineLineRelation::Skew);
  const Vector3 gap = r.pointOnB - r.pointOnA;
  REQUIRE_NEAR(dot(gap, a.direction()), 0.0, 1e-13);
  REQUIRE_NEAR(dot(gap, b.direction()), 0.0, 1e-13);
  REQUIRE_NEAR(norm(gap), r.distance, 1e-13);
  // The reported pair lies on the lines.
  REQUIRE_POINT_NEAR(a.pointAt(r.parameterA), r.pointOnA, 1e-14);
  REQUIRE_POINT_NEAR(b.pointAt(r.parameterB), r.pointOnB, 1e-14);
  // No sampled pair of points is closer than the reported distance.
  for (double s = -5; s <= 5; s += 0.5) {
    for (double t = -5; t <= 5; t += 0.5) {
      REQUIRE(distance(a.pointAt(s), b.pointAt(t)) >= r.distance - 1e-12);
    }
  }
}

TEST_CASE("Closest points: intersecting, parallel, coincident", "[geometry][lines]") {
  const Line3 x = Line3::fromPointDirection({0, 0, 0}, {1, 0, 0});

  const Line3 crossing = Line3::fromPointDirection({3, -5, 0}, {0, 1, 0});
  const auto hit = closestPointsBetweenLines(x, crossing);
  REQUIRE(hit.relation == LineLineRelation::Intersecting);
  REQUIRE_NEAR(hit.distance, 0.0, 1e-12);
  REQUIRE_POINT_NEAR(hit.pointOnA, (Point3{3, 0, 0}), 1e-13);
  REQUIRE_POINT_NEAR(hit.pointOnB, (Point3{3, 0, 0}), 1e-13);

  const Line3 parallel = Line3::fromPointDirection({7, 3, 4}, {-2, 0, 0});  // antiparallel dir
  const auto par = closestPointsBetweenLines(x, parallel);
  REQUIRE(par.relation == LineLineRelation::Parallel);
  REQUIRE_NEAR(par.distance, 5.0, 1e-14);
  REQUIRE_NEAR(norm(par.pointOnB - par.pointOnA), 5.0, 1e-14);  // representative pair is closest

  const Line3 same = Line3::fromPointDirection({-9, 0, 0}, {3, 0, 0});
  const auto co = closestPointsBetweenLines(x, same);
  REQUIRE(co.relation == LineLineRelation::Coincident);
  REQUIRE_NEAR(co.distance, 0.0, 1e-14);
}

TEST_CASE("Closest points: near-parallel lines are classified by the angular tolerance",
          "[geometry][lines][numerical]") {
  const Line3 a = Line3::fromPointDirection({0, 0, 0}, {1, 0, 0});
  // Offset by 1 in z; direction tilted by 1e-6 rad in the xy plane -> skew, distance 1.
  const double eps = 1e-6;
  const Line3 b = Line3::fromPointDirection({1000, 1, 1}, {std::cos(eps), std::sin(eps), 0});
  const auto skew = closestPointsBetweenLines(a, b);
  REQUIRE(skew.relation == LineLineRelation::Skew);
  REQUIRE_NEAR(skew.distance, 1.0, 1e-8);
  REQUIRE_NEAR(skew.pointOnA.z, 0.0, 1e-8);
  REQUIRE_NEAR(skew.pointOnB.z, 1.0, 1e-8);
  // The closest point is ~1e6 away from the origins: well-defined but far.
  REQUIRE_THAT(std::abs(skew.pointOnA.x), Catch::Matchers::WithinRel(1.0 / std::tan(eps) - 1000.0,
                                                                     1e-6));

  // Tilted by 1e-12 rad: below the angular tolerance -> reported parallel, distance still 1.
  const double tiny = 1e-12;
  const Line3 c = Line3::fromPointDirection({1000, 1, 1}, {std::cos(tiny), std::sin(tiny), 0});
  const auto par = closestPointsBetweenLines(a, c);
  REQUIRE(par.relation == LineLineRelation::Parallel);
  REQUIRE_THAT(par.distance, Catch::Matchers::WithinRel(std::sqrt(2.0), 1e-9));  // |(1,1)| offset
}

TEST_CASE("Geometry is covariant under rigid transforms", "[geometry][transform]") {
  const Line3 a = Line3::fromPointDirection({1, -2, 0.5}, {1, 2, 3});
  const Line3 b = Line3::fromPointDirection({-4, 3, 2}, {-2, 1, 1});
  const RigidTransform3 t = RigidTransform3::rotationAbout({2, 2, 2}, {1, -1, 2}, 0.9) *
                            RigidTransform3::translation({10, -3, 4});
  const Line3 ta = Line3::fromPointDirection(t.apply(a.origin()), t.apply(a.direction()));
  const Line3 tb = Line3::fromPointDirection(t.apply(b.origin()), t.apply(b.direction()));
  REQUIRE_NEAR(closestPointsBetweenLines(ta, tb).distance, closestPointsBetweenLines(a, b).distance,
               1e-12);

  const Plane3 plane = Plane3::fromPointNormal({0, 1, 0}, {1, 1, 1});
  const Plane3 tplane = Plane3::fromPointNormal(t.apply(plane.origin()), t.apply(plane.normal()));
  const Point3 p{3, -1, 6};
  REQUIRE_NEAR(distancePointToPlane(t.apply(p), tplane), distancePointToPlane(p, plane), 1e-12);
  const auto h = intersectLinePlane(a, plane);
  const auto th = intersectLinePlane(ta, tplane);
  REQUIRE(h.relation == LinePlaneRelation::UniquePoint);
  REQUIRE(th.relation == LinePlaneRelation::UniquePoint);
  REQUIRE_POINT_NEAR(t.apply(*h.point), *th.point, 1e-12);
}
